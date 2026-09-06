"""Ring-topology delta audit for forward validation.

The audit is intentionally graph-level and deterministic. It does not try to
prove a synthesis. It reports whether a proposed one-step transform creates,
removes, merges, or rewrites ring systems in ways that require stronger
mechanistic topology evidence before commit.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS

from ._rdkit_utils import canonical, parse_mol
from .mol_info import analyze_molecule
from .validation_findings import make_finding


def _ring_systems(atom_rings: Sequence[Set[int]]) -> List[Set[int]]:
    systems: List[Set[int]] = []
    for ring_atoms in atom_rings:
        merged = set(ring_atoms)
        keep: List[Set[int]] = []
        for system in systems:
            if merged & system:
                merged |= system
            else:
                keep.append(system)
        keep.append(merged)
        systems = keep

    changed = True
    while changed:
        changed = False
        compact: List[Set[int]] = []
        for system in systems:
            for idx, existing in enumerate(compact):
                if system & existing:
                    compact[idx] = existing | system
                    changed = True
                    break
            else:
                compact.append(set(system))
        systems = compact
    return systems


def _ring_features(smiles: str) -> Dict[str, Any]:
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    ri = mol.GetRingInfo()
    atom_rings = [set(ring) for ring in ri.AtomRings()]
    ring_sizes = sorted(len(ring) for ring in atom_rings)
    systems = _ring_systems(atom_rings)

    fused_pairs = 0
    bridged_pairs = 0
    spiro_pairs = 0
    for i in range(len(atom_rings)):
        for j in range(i + 1, len(atom_rings)):
            shared = len(atom_rings[i] & atom_rings[j])
            if shared >= 3:
                bridged_pairs += 1
            elif shared == 2:
                fused_pairs += 1
            elif shared == 1:
                spiro_pairs += 1

    junction_atoms = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if ri.NumAtomRings(atom.GetIdx()) > 1
    ]
    scaffold = analyze_molecule(smiles).get("scaffold", {})

    return {
        "ok": True,
        "smiles": canonical(smiles) or smiles,
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "ring_count": len(atom_rings),
        "ring_sizes": ring_sizes,
        "ring_system_sizes": sorted(len(system) for system in systems),
        "largest_ring_system_size": max((len(system) for system in systems), default=0),
        "ring_bond_count": sum(1 for bond in mol.GetBonds() if bond.IsInRing()),
        "fused_pairs": fused_pairs,
        "bridged_pairs": bridged_pairs,
        "spiro_pairs": spiro_pairs,
        "junction_atom_count": len(junction_atoms),
        "junction_atoms": junction_atoms,
        "tags": list(scaffold.get("tags", []) or []),
    }


def _mcs_coverage(product_mol: Chem.Mol, precursor_mol: Chem.Mol) -> float:
    try:
        mcs = rdFMCS.FindMCS(
            [product_mol, precursor_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
            timeout=2,
        )
        if not getattr(mcs, "smartsString", ""):
            return 0.0
        query = Chem.MolFromSmarts(mcs.smartsString)
        if query is None:
            return 0.0
        return round(query.GetNumHeavyAtoms() / max(product_mol.GetNumHeavyAtoms(), 1), 4)
    except Exception:
        return 0.0


def _select_major_precursor(product_smiles: str, precursors: Sequence[str]) -> Dict[str, Any]:
    product_mol = parse_mol(product_smiles)
    if product_mol is None:
        return {"index": None, "smiles": "", "mcs_coverage": 0.0}

    best: Optional[Dict[str, Any]] = None
    for idx, smi in enumerate(precursors):
        mol = parse_mol(smi)
        if mol is None:
            continue
        coverage = _mcs_coverage(product_mol, mol)
        features = _ring_features(smi)
        score = coverage * 100.0 + mol.GetNumHeavyAtoms()
        score += float(features.get("ring_count", 0)) * 5.0
        row = {
            "index": idx,
            "smiles": canonical(smi) or smi,
            "mcs_coverage": coverage,
            "score": score,
            "features": features,
        }
        if best is None or row["score"] > best["score"]:
            best = row

    if best is None:
        return {"index": None, "smiles": "", "mcs_coverage": 0.0}
    return best


def _combined_precursor_max(precursor_features: Sequence[Dict[str, Any]], key: str) -> int:
    values = [
        int(features.get(key, 0) or 0)
        for features in precursor_features
        if features.get("ok", True)
    ]
    return max(values, default=0)


def _combined_precursor_sum(precursor_features: Sequence[Dict[str, Any]], key: str) -> int:
    return sum(
        int(features.get(key, 0) or 0)
        for features in precursor_features
        if features.get("ok", True)
    )


def _combined_ring_sizes(precursor_features: Sequence[Dict[str, Any]]) -> Set[int]:
    sizes: Set[int] = set()
    for features in precursor_features:
        for size in features.get("ring_sizes", []) or []:
            try:
                sizes.add(int(size))
            except (TypeError, ValueError):
                continue
    return sizes


def _risk_level(delta_codes: Sequence[str]) -> str:
    critical = {
        "new_fused_ring_system",
        "new_spiro_center",
        "new_bridged_system",
        "new_macrocycle",
        "new_medium_ring",
        "ring_system_merge",
    }
    if any(code in critical for code in delta_codes):
        return "critical"
    if "ring_count_increase" in delta_codes or "ring_count_decrease" in delta_codes:
        return "high"
    if delta_codes:
        return "medium"
    return "none"


def _violation_from_delta(delta_codes: Sequence[str]) -> List[Dict[str, Any]]:
    violations: List[Dict[str, Any]] = []
    if "new_fused_ring_system" in delta_codes and "new_medium_ring" in delta_codes:
        violations.append(make_finding(
            code="new_fused_medium_ring_requires_evidence",
            severity="requires_evidence",
            source="ring_topology",
            message="a new fused medium-ring topology was detected",
            evidence={"delta_codes": list(delta_codes)},
            required_evidence=[
                "atom-mapped ring-atom source proof",
                "mechanistic topology proof",
            ],
        ))
    elif "new_fused_ring_system" in delta_codes:
        violations.append(make_finding(
            code="new_fused_ring_requires_evidence",
            severity="requires_evidence",
            source="ring_topology",
            message="a new fused-ring topology was detected",
            evidence={"delta_codes": list(delta_codes)},
            required_evidence=[
                "atom-mapped ring-atom source proof",
                "mechanistic topology proof",
            ],
        ))
    if "new_spiro_center" in delta_codes:
        violations.append(make_finding(
            code="new_spiro_center_requires_evidence",
            severity="requires_evidence",
            source="ring_topology",
            message="a new spiro junction was detected",
            evidence={"delta_codes": list(delta_codes)},
            required_evidence=[
                "atom-mapped spiro-center source proof",
                "mechanistic topology proof",
            ],
        ))
    if "new_bridged_system" in delta_codes:
        violations.append(make_finding(
            code="new_bridged_system_requires_evidence",
            severity="requires_evidence",
            source="ring_topology",
            message="a new bridged-ring topology was detected",
            evidence={"delta_codes": list(delta_codes)},
            required_evidence=[
                "atom-mapped bridgehead source proof",
                "mechanistic topology proof",
            ],
        ))
    if "new_macrocycle" in delta_codes:
        violations.append(make_finding(
            code="new_macrocycle_requires_evidence",
            severity="requires_evidence",
            source="ring_topology",
            message="a new macrocycle was detected",
            evidence={"delta_codes": list(delta_codes)},
            required_evidence=[
                "intramolecular tether proof",
                "mechanistic macrocyclization proof",
            ],
        ))
    return violations


def audit_ring_topology(
    product_smiles: str,
    precursor_smiles: Sequence[str],
    *,
    reaction_category: Optional[str] = None,
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compare product and precursor ring topology.

    The returned payload is gate-friendly: compact delta codes plus small
    feature summaries. It deliberately avoids large atom maps.
    """
    product_can = canonical(product_smiles) or product_smiles
    product = _ring_features(product_can)
    if not product.get("ok", True):
        return product

    precursor_features = [
        _ring_features(canonical(smi) or smi)
        for smi in precursor_smiles
        if parse_mol(smi) is not None
    ]
    major = _select_major_precursor(product_can, precursor_smiles)
    major_features = major.get("features") or {}

    inherited_prec_ring_count = _combined_precursor_sum(precursor_features, "ring_count")
    max_prec_fused = _combined_precursor_max(precursor_features, "fused_pairs")
    max_prec_spiro = _combined_precursor_max(precursor_features, "spiro_pairs")
    max_prec_bridged = _combined_precursor_max(precursor_features, "bridged_pairs")
    max_prec_junction = _combined_precursor_max(precursor_features, "junction_atom_count")
    max_prec_system = _combined_precursor_max(precursor_features, "largest_ring_system_size")
    precursor_ring_sizes = _combined_ring_sizes(precursor_features)

    delta_codes: List[str] = []
    if int(product.get("ring_count", 0)) > inherited_prec_ring_count:
        delta_codes.append("ring_count_increase")
    if int(product.get("ring_count", 0)) < inherited_prec_ring_count:
        delta_codes.append("ring_count_decrease")
    if int(product.get("fused_pairs", 0)) > max_prec_fused:
        delta_codes.append("new_fused_ring_system")
    if int(product.get("spiro_pairs", 0)) > max_prec_spiro:
        delta_codes.append("new_spiro_center")
    if int(product.get("bridged_pairs", 0)) > max_prec_bridged:
        delta_codes.append("new_bridged_system")
    if int(product.get("junction_atom_count", 0)) > max_prec_junction:
        delta_codes.append("new_junction_atoms")
    if int(product.get("largest_ring_system_size", 0)) > max_prec_system:
        delta_codes.append("ring_system_merge")
    if any(int(size) >= 12 for size in product.get("ring_sizes", []) or []) and not any(
        int(size) >= 12 for size in precursor_ring_sizes
    ):
        delta_codes.append("new_macrocycle")
    if any(7 <= int(size) <= 11 for size in product.get("ring_sizes", []) or []) and not any(
        7 <= int(size) <= 11 for size in precursor_ring_sizes
    ):
        delta_codes.append("new_medium_ring")
    if any(int(size) <= 4 for size in product.get("ring_sizes", []) or []) and not any(
        int(size) <= 4 for size in precursor_ring_sizes
    ):
        delta_codes.append("new_strained_small_ring")

    action_context = dict(action_context or {})
    if action_context.get("in_ring") and "ring_bond_edit" not in delta_codes:
        delta_codes.append("ring_bond_edit")

    # Preserve order while removing duplicates.
    delta_codes = list(dict.fromkeys(delta_codes))
    risk_level = _risk_level(delta_codes)
    violations = _violation_from_delta(delta_codes)

    applies = bool(
        product.get("ring_count")
        or precursor_features
        or action_context.get("in_ring")
        or delta_codes
    )
    return {
        "ok": True,
        "applies": applies,
        "risk_level": risk_level,
        "reaction_category": reaction_category or "",
        "major_precursor_index": major.get("index"),
        "major_precursor": major.get("smiles", ""),
        "major_precursor_mcs_coverage": major.get("mcs_coverage", 0.0),
        "product": product,
        "major_precursor_features": major_features,
        "precursor_max": {
            "ring_count": inherited_prec_ring_count,
            "fused_pairs": max_prec_fused,
            "spiro_pairs": max_prec_spiro,
            "bridged_pairs": max_prec_bridged,
            "junction_atom_count": max_prec_junction,
            "largest_ring_system_size": max_prec_system,
            "ring_sizes": sorted(precursor_ring_sizes),
        },
        "delta_codes": delta_codes,
        "violations": violations,
        "required_evidence": [
            "atom-mapped ring-atom source proof",
            "mechanistic topology proof",
        ] if risk_level in {"high", "critical"} else [],
        "summary": (
            "no ring topology delta"
            if not delta_codes
            else "ring topology delta: " + ", ".join(delta_codes)
        ),
    }
