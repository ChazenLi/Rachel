"""Generic graph-delta audit for forward validation.

This module reports graph facts: fragment merging, mapped heavy-bond changes,
topological distance changes, and low MCS coverage. It deliberately does not
decide whether a reaction should commit.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS

from ._rdkit_utils import canonical, parse_mol
from .validation_findings import make_finding


def _combine_mols(mols: Sequence[Chem.Mol]) -> Optional[Chem.Mol]:
    if not mols:
        return None
    combined = mols[0]
    for mol in mols[1:]:
        combined = Chem.CombineMols(combined, mol)
    return combined


def _frag_index_by_atom(mol: Chem.Mol) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for frag_idx, atom_ids in enumerate(Chem.GetMolFrags(mol)):
        for atom_idx in atom_ids:
            out[int(atom_idx)] = frag_idx
    return out


def _bond_evidence(mol: Chem.Mol, atom_a: int, atom_b: int) -> Optional[Dict[str, Any]]:
    bond = mol.GetBondBetweenAtoms(int(atom_a), int(atom_b))
    if bond is None:
        return None
    return {
        "order": float(bond.GetBondTypeAsDouble()),
        "aromatic": bool(bond.GetIsAromatic()),
        "in_ring": bool(bond.IsInRing()),
    }


def _sample_append(rows: List[Dict[str, Any]], row: Dict[str, Any], limit: int = 8) -> None:
    if len(rows) < limit:
        rows.append(row)


def _risk_level(delta_codes: Sequence[str]) -> str:
    if "apparent_scaffold_jump" in delta_codes:
        return "high"
    medium = {
        "mapped_heavy_bond_formation",
        "mapped_heavy_bond_cleavage",
        "mapped_topological_distance_change",
    }
    if any(code in medium for code in delta_codes):
        return "medium"
    if delta_codes:
        return "low"
    return "none"


def audit_graph_delta(
    product_smiles: str,
    precursor_smiles: Sequence[str],
    *,
    reaction_category: Optional[str] = None,
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compare product graph to the combined precursor graph.

    The output is fact-only. Later policy code decides how to interpret these
    deltas for a particular reaction family.
    """
    product = parse_mol(product_smiles)
    precursor_mols = [parse_mol(smi) for smi in precursor_smiles]
    precursor_mols = [mol for mol in precursor_mols if mol is not None]
    if product is None:
        return {"ok": False, "error": "invalid product SMILES", "input": product_smiles}
    combined = _combine_mols(precursor_mols)
    if combined is None:
        return {"ok": False, "error": "no valid precursor molecules"}

    product_heavy = product.GetNumHeavyAtoms()
    precursor_heavy = combined.GetNumHeavyAtoms()
    product_frag_count = len(Chem.GetMolFrags(product))
    precursor_frag_count = sum(len(Chem.GetMolFrags(mol)) for mol in precursor_mols)

    try:
        mcs = rdFMCS.FindMCS(
            [product, combined],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=3,
        )
    except Exception as exc:
        return {
            "ok": True,
            "applies": True,
            "reaction_category": reaction_category or "",
            "delta_codes": ["mcs_failed"],
            "risk_level": "low",
            "findings": [
                make_finding(
                    code="mcs_failed",
                    severity="warning",
                    source="graph_delta",
                    message="MCS graph-delta comparison failed",
                    evidence={"error": str(exc)},
                )
            ],
            "violations": [],
            "summary": "graph delta MCS failed",
        }

    mcs_num_atoms = int(getattr(mcs, "numAtoms", 0) or 0)
    mcs_smarts = getattr(mcs, "smartsString", "") or ""
    mcs_product_coverage = round(mcs_num_atoms / max(product_heavy, 1), 4)
    mcs_precursor_coverage = round(mcs_num_atoms / max(precursor_heavy, 1), 4)

    delta_codes: List[str] = []
    findings: List[Dict[str, Any]] = []

    if precursor_frag_count > 1 and product_frag_count == 1:
        delta_codes.append("fragment_merge")
        findings.append(make_finding(
            code="fragment_merge",
            severity="info",
            source="graph_delta",
            message="multiple precursor fragments map to a single product component",
            evidence={
                "precursor_fragment_count": precursor_frag_count,
                "product_fragment_count": product_frag_count,
            },
        ))

    if mcs_product_coverage < 0.35 and product_heavy >= 8:
        delta_codes.append("apparent_scaffold_jump")
        findings.append(make_finding(
            code="apparent_scaffold_jump",
            severity="requires_evidence",
            source="graph_delta",
            message="MCS covers little of the product graph; this may indicate scaffold jump or multi-step compression",
            evidence={
                "mcs_num_atoms": mcs_num_atoms,
                "mcs_product_coverage": mcs_product_coverage,
                "product_heavy_atoms": product_heavy,
            },
            required_evidence=[
                "mechanistic scaffold-editing proof",
                "atom-mapped source of retained scaffold atoms",
            ],
        ))
    elif mcs_product_coverage < 0.5 and product_heavy >= 8:
        delta_codes.append("low_mcs_coverage")
        findings.append(make_finding(
            code="low_mcs_coverage",
            severity="warning",
            source="graph_delta",
            message="MCS coverage of product graph is low",
            evidence={
                "mcs_num_atoms": mcs_num_atoms,
                "mcs_product_coverage": mcs_product_coverage,
                "product_heavy_atoms": product_heavy,
            },
        ))

    query = Chem.MolFromSmarts(mcs_smarts) if mcs_smarts else None
    match_product: Tuple[int, ...] = tuple(product.GetSubstructMatch(query)) if query is not None else tuple()
    match_precursor: Tuple[int, ...] = tuple(combined.GetSubstructMatch(query)) if query is not None else tuple()

    new_bonds: List[Dict[str, Any]] = []
    broken_bonds: List[Dict[str, Any]] = []
    order_changes: List[Dict[str, Any]] = []
    distance_changes: List[Dict[str, Any]] = []

    if query is not None and len(match_product) == len(match_precursor) and match_product:
        p_to_r = {int(p): int(r) for p, r in zip(match_product, match_precursor)}
        r_to_p = {int(r): int(p) for p, r in zip(match_product, match_precursor)}

        for bond in product.GetBonds():
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()
            if a not in p_to_r or b not in p_to_r:
                continue
            r_a = p_to_r[a]
            r_b = p_to_r[b]
            precursor_bond = _bond_evidence(combined, r_a, r_b)
            product_bond = _bond_evidence(product, a, b)
            if precursor_bond is None:
                _sample_append(new_bonds, {
                    "product_atoms": [int(a), int(b)],
                    "precursor_atoms": [int(r_a), int(r_b)],
                    "product_bond": product_bond,
                })
            elif product_bond and product_bond != precursor_bond:
                _sample_append(order_changes, {
                    "product_atoms": [int(a), int(b)],
                    "precursor_atoms": [int(r_a), int(r_b)],
                    "product_bond": product_bond,
                    "precursor_bond": precursor_bond,
                })

        for bond in combined.GetBonds():
            r_a = bond.GetBeginAtomIdx()
            r_b = bond.GetEndAtomIdx()
            if r_a not in r_to_p or r_b not in r_to_p:
                continue
            p_a = r_to_p[r_a]
            p_b = r_to_p[r_b]
            if product.GetBondBetweenAtoms(p_a, p_b) is None:
                _sample_append(broken_bonds, {
                    "precursor_atoms": [int(r_a), int(r_b)],
                    "product_atoms": [int(p_a), int(p_b)],
                    "precursor_bond": _bond_evidence(combined, r_a, r_b),
                })

        precursor_frag_by_atom = _frag_index_by_atom(combined)
        dist_product = Chem.GetDistanceMatrix(product)
        dist_precursor = Chem.GetDistanceMatrix(combined)
        mapped_product = list(p_to_r.keys())
        for idx, p_a in enumerate(mapped_product):
            for p_b in mapped_product[idx + 1:]:
                r_a = p_to_r[p_a]
                r_b = p_to_r[p_b]
                if precursor_frag_by_atom.get(r_a) != precursor_frag_by_atom.get(r_b):
                    continue
                d_product = int(dist_product[p_a][p_b])
                d_precursor = int(dist_precursor[r_a][r_b])
                if d_product != d_precursor:
                    _sample_append(distance_changes, {
                        "product_atoms": [int(p_a), int(p_b)],
                        "precursor_atoms": [int(r_a), int(r_b)],
                        "product_distance": d_product,
                        "precursor_distance": d_precursor,
                    })

    if new_bonds:
        delta_codes.append("mapped_heavy_bond_formation")
        findings.append(make_finding(
            code="mapped_heavy_bond_formation",
            severity="warning",
            source="graph_delta",
            message="mapped heavy atoms are bonded in the product but not in the precursor graph",
            evidence={"sample": new_bonds, "count_sampled": len(new_bonds)},
        ))
    if broken_bonds:
        delta_codes.append("mapped_heavy_bond_cleavage")
        findings.append(make_finding(
            code="mapped_heavy_bond_cleavage",
            severity="warning",
            source="graph_delta",
            message="mapped heavy-atom bonds from precursors are absent in the product graph",
            evidence={"sample": broken_bonds, "count_sampled": len(broken_bonds)},
        ))
    if order_changes:
        delta_codes.append("mapped_bond_order_change")
        findings.append(make_finding(
            code="mapped_bond_order_change",
            severity="info",
            source="graph_delta",
            message="mapped heavy-atom bond attributes changed",
            evidence={"sample": order_changes, "count_sampled": len(order_changes)},
        ))
    if distance_changes:
        delta_codes.append("mapped_topological_distance_change")
        findings.append(make_finding(
            code="mapped_topological_distance_change",
            severity="warning",
            source="graph_delta",
            message="topological distances between retained mapped atoms changed",
            evidence={"sample": distance_changes, "count_sampled": len(distance_changes)},
            required_evidence=[
                "intended ring opening/closure/rearrangement explanation",
                "mechanistic topology-distance proof",
            ],
        ))

    delta_codes = list(dict.fromkeys(delta_codes))
    risk_level = _risk_level(delta_codes)
    action_context = dict(action_context or {})

    return {
        "ok": True,
        "applies": bool(delta_codes or mcs_num_atoms),
        "source": "graph_delta",
        "risk_level": risk_level,
        "reaction_category": reaction_category or "",
        "product": canonical(product_smiles) or product_smiles,
        "precursors": [canonical(smi) or smi for smi in precursor_smiles],
        "product_heavy_atoms": product_heavy,
        "combined_precursor_heavy_atoms": precursor_heavy,
        "product_fragment_count": product_frag_count,
        "precursor_fragment_count": precursor_frag_count,
        "mcs_num_atoms": mcs_num_atoms,
        "mcs_product_coverage": mcs_product_coverage,
        "mcs_precursor_coverage": mcs_precursor_coverage,
        "delta_codes": delta_codes,
        "findings": findings,
        "violations": [
            finding for finding in findings
            if finding.get("severity") in {"requires_evidence", "hard_fail"}
        ],
        "action_context": action_context,
        "summary": (
            "no graph delta"
            if not delta_codes
            else "graph delta: " + ", ".join(delta_codes)
        ),
    }
