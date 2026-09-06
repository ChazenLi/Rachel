"""MCS-derived atom-source audit for forward validation.

This is not a true experimentally atom-mapped reaction. It is a deterministic
RDKit proof aid: build a best-effort product-to-precursor atom correspondence,
then report preserved atoms, mapped bond formations/cleavages, ambiguity, and
missing proof fields for LLM review.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS

from ._rdkit_utils import canonical, parse_mol
from .validation_findings import make_finding
from .topology_intent import action_declares_topology_change


def _combine_mols(mols: Sequence[Chem.Mol]) -> Optional[Chem.Mol]:
    if not mols:
        return None
    combined = mols[0]
    for mol in mols[1:]:
        combined = Chem.CombineMols(combined, mol)
    return combined


def _precursor_atom_index(mols: Sequence[Chem.Mol]) -> Dict[int, Dict[str, int]]:
    out: Dict[int, Dict[str, int]] = {}
    offset = 0
    for precursor_idx, mol in enumerate(mols):
        for atom in mol.GetAtoms():
            out[offset + atom.GetIdx()] = {
                "precursor_index": int(precursor_idx),
                "precursor_atom": int(atom.GetIdx()),
            }
        offset += mol.GetNumAtoms()
    return out


def _frag_index_by_atom(mol: Chem.Mol) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for frag_idx, atom_ids in enumerate(Chem.GetMolFrags(mol)):
        for atom_idx in atom_ids:
            out[int(atom_idx)] = frag_idx
    return out


def _bond_row(
    product: Chem.Mol,
    combined: Chem.Mol,
    p_to_r: Mapping[int, int],
    precursor_lookup: Mapping[int, Dict[str, int]],
    precursor_frag_by_atom: Mapping[int, int],
    atom_a: int,
    atom_b: int,
) -> Dict[str, Any]:
    r_a = int(p_to_r[int(atom_a)])
    r_b = int(p_to_r[int(atom_b)])
    product_bond = product.GetBondBetweenAtoms(int(atom_a), int(atom_b))
    precursor_bond = combined.GetBondBetweenAtoms(r_a, r_b)
    precursor_a = precursor_lookup.get(r_a, {})
    precursor_b = precursor_lookup.get(r_b, {})
    row: Dict[str, Any] = {
        "product_atoms": [int(atom_a), int(atom_b)],
        "product_symbols": [
            product.GetAtomWithIdx(int(atom_a)).GetSymbol(),
            product.GetAtomWithIdx(int(atom_b)).GetSymbol(),
        ],
        "precursor_atoms": [r_a, r_b],
        "precursor_indices": [
            precursor_a.get("precursor_index"),
            precursor_b.get("precursor_index"),
        ],
        "precursor_local_atoms": [
            precursor_a.get("precursor_atom"),
            precursor_b.get("precursor_atom"),
        ],
        "precursor_fragments": [
            precursor_frag_by_atom.get(r_a),
            precursor_frag_by_atom.get(r_b),
        ],
    }
    if product_bond is not None:
        row["product_bond_order"] = float(product_bond.GetBondTypeAsDouble())
        row["product_bond_in_ring"] = bool(product_bond.IsInRing())
    if precursor_bond is not None:
        row["precursor_bond_order"] = float(precursor_bond.GetBondTypeAsDouble())
        row["precursor_bond_in_ring"] = bool(precursor_bond.IsInRing())
    row["same_precursor"] = row["precursor_indices"][0] == row["precursor_indices"][1]
    row["same_precursor_fragment"] = row["precursor_fragments"][0] == row["precursor_fragments"][1]
    return row


def _sample(rows: List[Dict[str, Any]], row: Dict[str, Any], limit: int = 8) -> None:
    if len(rows) < limit:
        rows.append(row)


def _count_unmapped(mol: Chem.Mol, mapped_atoms: set[int]) -> Dict[str, int]:
    counter = Counter(
        atom.GetSymbol()
        for atom in mol.GetAtoms()
        if atom.GetIdx() not in mapped_atoms and atom.GetAtomicNum() > 1
    )
    return dict(sorted(counter.items()))


def _anchor_candidates(product: Chem.Mol, p_to_r: Mapping[int, int], limit: int = 12) -> List[Dict[str, Any]]:
    ri = product.GetRingInfo()
    rows: List[Dict[str, Any]] = []
    for p_idx, r_idx in p_to_r.items():
        atom = product.GetAtomWithIdx(int(p_idx))
        symbol = atom.GetSymbol()
        ring_count = int(ri.NumAtomRings(int(p_idx)))
        if symbol == "C" and ring_count == 0 and atom.GetDegree() < 3:
            continue
        rows.append({
            "product_atom": int(p_idx),
            "precursor_atom": int(r_idx),
            "symbol": symbol,
            "ring_count": ring_count,
            "degree": int(atom.GetDegree()),
        })
        if len(rows) >= limit:
            break
    return rows


def _custom_topology_context(action_context: Mapping[str, Any]) -> bool:
    source = str(action_context.get("source", "") or "")
    if source not in {"custom_precursors", "llm_proposed"}:
        return False
    return bool(action_context.get("in_ring")) or action_declares_topology_change(action_context)


def _missing_declared_fields(action_context: Mapping[str, Any]) -> List[str]:
    required = [
        "intended_deltas",
        "expected_ring_change",
        "changed_bonds",
        "preserved_anchors",
        "mechanistic_evidence",
    ]
    return [key for key in required if action_context.get(key) in (None, "", [], {})]


def audit_atom_mapping(
    product_smiles: str,
    precursor_smiles: Sequence[str],
    *,
    reaction_category: Optional[str] = None,
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return compact MCS-derived atom-source evidence."""
    action_context = dict(action_context or {})
    product = parse_mol(product_smiles)
    precursor_mols = [parse_mol(smi) for smi in precursor_smiles]
    precursor_mols = [mol for mol in precursor_mols if mol is not None]
    if product is None:
        return {"ok": False, "error": "invalid product SMILES", "input": product_smiles}
    combined = _combine_mols(precursor_mols)
    if combined is None:
        return {"ok": False, "error": "no valid precursor molecules"}

    findings: List[Dict[str, Any]] = []
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
            "source": "atom_mapping",
            "status": "mcs_failed",
            "confidence": "low",
            "reaction_category": reaction_category or "",
            "delta_codes": ["atom_map_failed"],
            "findings": [
                make_finding(
                    code="atom_map_failed",
                    severity="warning",
                    source="atom_mapping",
                    message="MCS atom-source audit failed",
                    evidence={"error": str(exc)},
                )
            ],
            "violations": [],
            "summary": "atom mapping failed",
        }

    mcs_num_atoms = int(getattr(mcs, "numAtoms", 0) or 0)
    product_heavy = product.GetNumHeavyAtoms()
    precursor_heavy = combined.GetNumHeavyAtoms()
    product_coverage = round(mcs_num_atoms / max(product_heavy, 1), 4)
    precursor_coverage = round(mcs_num_atoms / max(precursor_heavy, 1), 4)
    smarts = getattr(mcs, "smartsString", "") or ""
    query = Chem.MolFromSmarts(smarts) if smarts else None
    if query is None or mcs_num_atoms == 0:
        findings.append(make_finding(
            code="atom_map_unavailable",
            severity="missing_evidence",
            source="atom_mapping",
            message="no MCS atom-source map could be built",
            evidence={"mcs_num_atoms": mcs_num_atoms},
            required_evidence=["manual atom-source proof for changed bonds"],
        ))
        return {
            "ok": True,
            "applies": True,
            "source": "atom_mapping",
            "status": "no_map",
            "confidence": "low",
            "reaction_category": reaction_category or "",
            "mcs_num_atoms": mcs_num_atoms,
            "mcs_product_coverage": product_coverage,
            "mcs_precursor_coverage": precursor_coverage,
            "delta_codes": ["atom_map_unavailable"],
            "findings": findings,
            "violations": findings,
            "summary": "no atom-source map",
        }

    product_matches = product.GetSubstructMatches(query, uniquify=True, maxMatches=16)
    precursor_matches = combined.GetSubstructMatches(query, uniquify=True, maxMatches=16)
    product_match = tuple(product_matches[0]) if product_matches else tuple()
    precursor_match = tuple(precursor_matches[0]) if precursor_matches else tuple()
    if not product_match or not precursor_match or len(product_match) != len(precursor_match):
        findings.append(make_finding(
            code="atom_map_match_failed",
            severity="missing_evidence",
            source="atom_mapping",
            message="MCS exists but product/precursor substructure match failed",
            evidence={
                "product_match_count": len(product_matches),
                "precursor_match_count": len(precursor_matches),
            },
            required_evidence=["manual atom-source proof for changed bonds"],
        ))
        return {
            "ok": True,
            "applies": True,
            "source": "atom_mapping",
            "status": "match_failed",
            "confidence": "low",
            "reaction_category": reaction_category or "",
            "mcs_num_atoms": mcs_num_atoms,
            "mcs_product_coverage": product_coverage,
            "mcs_precursor_coverage": precursor_coverage,
            "delta_codes": ["atom_map_unavailable"],
            "findings": findings,
            "violations": findings,
            "summary": "atom-source map match failed",
        }

    p_to_r = {int(p): int(r) for p, r in zip(product_match, precursor_match)}
    r_to_p = {int(r): int(p) for p, r in p_to_r.items()}
    precursor_lookup = _precursor_atom_index(precursor_mols)
    precursor_frag_by_atom = _frag_index_by_atom(combined)

    formed_bonds: List[Dict[str, Any]] = []
    cleaved_bonds: List[Dict[str, Any]] = []
    bond_order_changes: List[Dict[str, Any]] = []

    for bond in product.GetBonds():
        p_a = int(bond.GetBeginAtomIdx())
        p_b = int(bond.GetEndAtomIdx())
        if p_a not in p_to_r or p_b not in p_to_r:
            continue
        r_a = p_to_r[p_a]
        r_b = p_to_r[p_b]
        precursor_bond = combined.GetBondBetweenAtoms(r_a, r_b)
        if precursor_bond is None:
            _sample(formed_bonds, _bond_row(product, combined, p_to_r, precursor_lookup, precursor_frag_by_atom, p_a, p_b))
        elif float(precursor_bond.GetBondTypeAsDouble()) != float(bond.GetBondTypeAsDouble()):
            _sample(bond_order_changes, _bond_row(product, combined, p_to_r, precursor_lookup, precursor_frag_by_atom, p_a, p_b))

    for bond in combined.GetBonds():
        r_a = int(bond.GetBeginAtomIdx())
        r_b = int(bond.GetEndAtomIdx())
        if r_a not in r_to_p or r_b not in r_to_p:
            continue
        p_a = r_to_p[r_a]
        p_b = r_to_p[r_b]
        if product.GetBondBetweenAtoms(p_a, p_b) is None:
            _sample(cleaved_bonds, _bond_row(product, combined, p_to_r, precursor_lookup, precursor_frag_by_atom, p_a, p_b))

    mapped_atoms = set(p_to_r)
    delta_codes: List[str] = []
    if formed_bonds:
        delta_codes.append("atom_mapped_bond_formation")
        findings.append(make_finding(
            code="atom_mapped_bond_formation",
            severity="info",
            source="atom_mapping",
            message="MCS-mapped atoms form new product bonds absent in the precursor graph",
            evidence={"sample": formed_bonds, "count_sampled": len(formed_bonds)},
        ))
    if cleaved_bonds:
        delta_codes.append("atom_mapped_bond_cleavage")
        findings.append(make_finding(
            code="atom_mapped_bond_cleavage",
            severity="info",
            source="atom_mapping",
            message="MCS-mapped precursor bonds are absent in the product graph",
            evidence={"sample": cleaved_bonds, "count_sampled": len(cleaved_bonds)},
        ))
    if bond_order_changes:
        delta_codes.append("atom_mapped_bond_order_change")
        findings.append(make_finding(
            code="atom_mapped_bond_order_change",
            severity="info",
            source="atom_mapping",
            message="MCS-mapped bonds changed order",
            evidence={"sample": bond_order_changes, "count_sampled": len(bond_order_changes)},
        ))

    ambiguous = len(product_matches) > 1 or len(precursor_matches) > 1
    if ambiguous:
        delta_codes.append("atom_map_ambiguous")
        findings.append(make_finding(
            code="atom_map_ambiguous",
            severity="warning",
            source="atom_mapping",
            message="multiple MCS matches exist; atom-source proof is ambiguous",
            evidence={
                "product_match_count": len(product_matches),
                "precursor_match_count": len(precursor_matches),
            },
            required_evidence=["manual disambiguation of changed-bond atom sources"],
        ))

    if product_coverage < 0.5 and product_heavy >= 8:
        delta_codes.append("atom_map_low_coverage")
        findings.append(make_finding(
            code="atom_map_low_coverage",
            severity="requires_evidence",
            source="atom_mapping",
            message="MCS atom map covers less than half of the product",
            evidence={
                "mcs_product_coverage": product_coverage,
                "mcs_num_atoms": mcs_num_atoms,
                "product_heavy_atoms": product_heavy,
            },
            required_evidence=["independent atom-source proof for retained scaffold atoms"],
        ))

    missing_fields = _missing_declared_fields(action_context) if _custom_topology_context(action_context) else []
    if missing_fields:
        delta_codes.append("atom_mapping_declared_evidence_missing")
        findings.append(make_finding(
            code="atom_mapping_declared_evidence_missing",
            severity="missing_evidence",
            source="atom_mapping",
            message="custom topology action lacks declared atom-source audit fields",
            evidence={"missing_fields": missing_fields},
            required_evidence=[
                "intended_deltas",
                "expected_ring_change",
                "changed_bonds",
                "preserved_anchors",
                "mechanistic_evidence",
            ],
        ))

    if product_coverage >= 0.75 and not ambiguous:
        confidence = "high"
    elif product_coverage >= 0.5:
        confidence = "medium"
    else:
        confidence = "low"
    if confidence == "high":
        status = "map_available"
    elif ambiguous:
        status = "ambiguous"
    elif confidence == "medium":
        status = "partial_map"
    else:
        status = "weak_map"

    mapped_sample = []
    for p_idx, r_idx in list(p_to_r.items())[:16]:
        precursor_info = precursor_lookup.get(r_idx, {})
        mapped_sample.append({
            "product_atom": int(p_idx),
            "product_symbol": product.GetAtomWithIdx(int(p_idx)).GetSymbol(),
            "precursor_atom": int(r_idx),
            "precursor_symbol": combined.GetAtomWithIdx(int(r_idx)).GetSymbol(),
            "precursor_index": precursor_info.get("precursor_index"),
            "precursor_local_atom": precursor_info.get("precursor_atom"),
            "precursor_fragment": precursor_frag_by_atom.get(int(r_idx)),
        })

    return {
        "ok": True,
        "applies": True,
        "source": "atom_mapping",
        "status": status,
        "confidence": confidence,
        "reaction_category": reaction_category or "",
        "product": canonical(product_smiles) or product_smiles,
        "precursors": [canonical(smi) or smi for smi in precursor_smiles],
        "mcs_num_atoms": mcs_num_atoms,
        "mcs_product_coverage": product_coverage,
        "mcs_precursor_coverage": precursor_coverage,
        "ambiguity": {
            "product_match_count": len(product_matches),
            "precursor_match_count": len(precursor_matches),
            "ambiguous": ambiguous,
        },
        "mapped_atom_count": len(p_to_r),
        "mapped_atom_sample": mapped_sample,
        "anchor_candidates": _anchor_candidates(product, p_to_r),
        "formed_bonds": formed_bonds,
        "cleaved_bonds": cleaved_bonds,
        "bond_order_changes": bond_order_changes,
        "unmapped_product_atoms": _count_unmapped(product, mapped_atoms),
        "unmapped_precursor_atoms": _count_unmapped(combined, set(r_to_p)),
        "declared_evidence": {
            "changed_bonds_provided": action_context.get("changed_bonds") not in (None, "", [], {}),
            "preserved_anchors_provided": action_context.get("preserved_anchors") not in (None, "", [], {}),
            "family_evidence_provided": action_context.get("family_evidence") not in (None, "", [], {}),
            "missing_fields": missing_fields,
        },
        "delta_codes": list(dict.fromkeys(delta_codes)),
        "findings": findings,
        "violations": [
            finding for finding in findings
            if finding.get("severity") in {"missing_evidence", "requires_evidence", "hard_fail"}
        ],
        "summary": (
            f"atom map {status}; coverage={product_coverage}; "
            f"formed={len(formed_bonds)} cleaved={len(cleaved_bonds)}"
        ),
    }
