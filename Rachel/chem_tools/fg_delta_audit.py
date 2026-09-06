"""Functional-group delta audit for forward validation.

This module reports functional-group changes as evidence. It is deliberately
not a reaction classifier and does not decide commit policy.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Optional, Sequence

from rdkit import Chem
from rdkit.Chem import rdFMCS

from ._rdkit_utils import canonical, parse_mol
from .fg_detect import detect_functional_groups, detect_protecting_groups
from .validation_findings import make_finding


_CARBONYL_GROUPS = {
    "carbonyl_generic",
    "aldehyde",
    "ketone",
    "cyclic_ketone",
}
_ALCOHOL_GROUPS = {
    "alcohol",
    "primary_alcohol",
    "secondary_alcohol",
    "tertiary_alcohol",
    "phenol",
}
_AMINE_GROUPS = {
    "primary_amine",
    "secondary_amine",
    "tertiary_amine",
    "aromatic_amine",
    "secondary_aromatic_amine",
}
_NITRO_GROUPS = {"nitro", "nitro_aromatic", "nitroalkane", "nitronate"}
_HALIDE_GROUPS = {
    "aryl_halide",
    "alkyl_halide",
    "alkyl_chloride",
    "alkyl_bromide",
    "alkyl_iodide",
    "alkyl_fluoride",
    "benzyl_halide",
    "allyl_halide",
    "vinyl_halide",
    "acyl_halide",
    "acyl_chloride",
}
_PROTECTION_GROUPS = {
    "boc_group",
    "cbz_group",
    "fmoc_group",
    "alloc_group",
    "benzyl_group",
    "pmb_group",
    "mom_ether",
    "mem_ether",
    "sem_ether",
    "thp_ether",
    "nosyl_group",
    "tosyl_on_N",
    "mesyl_on_N",
    "acetonide_on_diol",
    "benzylidene_acetal",
    "silyl_ether_tms",
    "silyl_ether_tes",
    "silyl_ether_tbs",
    "silyl_ether_tips",
}
_EXPLICIT_PG_NAMES = {
    "Boc",
    "Cbz",
    "Fmoc",
    "Alloc",
    "Troc",
    "Teoc",
    "Tosyl_N",
    "Nosyl_N",
    "Phthaloyl_N",
    "PMB_N",
    "Dde_N",
    "Mtt_N",
    "Bn_N",
    "TBS",
    "TBDPS",
    "TIPS",
    "TMS",
    "TES",
    "Bn_O",
    "PMB_O",
    "DMB_O",
    "Nap_O",
    "THP",
    "MOM",
    "MEM",
    "SEM",
    "Acetonide",
    "Benzylidene",
    "Trityl_O",
    "MMTr_O",
    "Mtt_S",
    "Trityl_S",
    "Acm_S",
    "StBu_S",
    "Mob_S",
    "Bn_S",
    "Acetal",
    "Ketal",
    "Dioxolane",
    "Dithiane",
    "Dithiolane",
}


def _group_counts(smiles: str) -> Counter:
    result = detect_functional_groups(smiles)
    if not result.get("ok"):
        return Counter()
    groups = result.get("groups", {}) or {}
    return Counter({
        name: int(info.get("count", 0) or 0)
        for name, info in groups.items()
        if isinstance(info, dict)
    })


def _protecting_group_counts(smiles: str) -> Counter:
    result = detect_protecting_groups(smiles)
    detected = result.get("detected", []) if isinstance(result, dict) else []
    return Counter({
        str(item.get("name", "")): int(item.get("count", 0) or 0)
        for item in detected
        if item.get("name") in _EXPLICIT_PG_NAMES
    })


def _halogen_atom_counts(smiles: str) -> Counter:
    mol = parse_mol(smiles)
    if mol is None:
        return Counter()
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() in {"F", "Cl", "Br", "I"})


def _mcs_coverage(product_mol: Chem.Mol, precursor_mol: Chem.Mol) -> float:
    try:
        mcs = rdFMCS.FindMCS(
            [product_mol, precursor_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=2,
        )
        return float(getattr(mcs, "numAtoms", 0) or 0) / max(product_mol.GetNumHeavyAtoms(), 1)
    except Exception:
        return 0.0


def _select_major_precursor(product_mol: Chem.Mol, precursor_smiles: Sequence[str]) -> Dict[str, Any]:
    best: Optional[Dict[str, Any]] = None
    for idx, smi in enumerate(precursor_smiles):
        mol = parse_mol(smi)
        if mol is None:
            continue
        coverage = _mcs_coverage(product_mol, mol)
        score = coverage * 100.0 + mol.GetNumHeavyAtoms()
        row = {
            "index": idx,
            "smiles": canonical(smi) or smi,
            "mcs_coverage": round(coverage, 4),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "score": score,
        }
        if best is None or row["score"] > best["score"]:
            best = row
    return best or {"index": None, "smiles": "", "mcs_coverage": 0.0, "heavy_atoms": 0}


def _sum_counts(smiles_list: Sequence[str], fn) -> Counter:
    total: Counter = Counter()
    for smi in smiles_list:
        total += fn(smi)
    return total


def _sum_for(counter: Counter, names: Sequence[str] | set[str]) -> int:
    return sum(int(counter.get(name, 0) or 0) for name in names)


def _delta_summary(product_counts: Counter, precursor_counts: Counter) -> Dict[str, Dict[str, int]]:
    keys = set(product_counts) | set(precursor_counts)
    increased: Dict[str, int] = {}
    decreased: Dict[str, int] = {}
    for key in sorted(keys):
        diff = int(product_counts.get(key, 0) or 0) - int(precursor_counts.get(key, 0) or 0)
        if diff > 0:
            increased[key] = diff
        elif diff < 0:
            decreased[key] = -diff
    return {"increased": increased, "decreased": decreased}


def audit_fg_delta(
    product_smiles: str,
    precursor_smiles: Sequence[str],
    *,
    reaction_category: Optional[str] = None,
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Compare product and precursor functional-group counts."""
    product_mol = parse_mol(product_smiles)
    if product_mol is None:
        return {"ok": False, "error": "invalid product SMILES", "input": product_smiles}
    valid_precursors = [smi for smi in precursor_smiles if parse_mol(smi) is not None]
    if not valid_precursors:
        return {"ok": False, "error": "no valid precursor molecules"}

    product_fg = _group_counts(product_smiles)
    combined_precursor_fg = _sum_counts(valid_precursors, _group_counts)
    product_pg = _protecting_group_counts(product_smiles)
    combined_precursor_pg = _sum_counts(valid_precursors, _protecting_group_counts)
    product_halogen_atoms = _halogen_atom_counts(product_smiles)
    combined_precursor_halogen_atoms = _sum_counts(valid_precursors, _halogen_atom_counts)

    major = _select_major_precursor(product_mol, valid_precursors)
    major_smiles = str(major.get("smiles") or "")
    major_fg = _group_counts(major_smiles) if major_smiles else Counter()
    major_pg = _protecting_group_counts(major_smiles) if major_smiles else Counter()
    major_halogen_atoms = _halogen_atom_counts(major_smiles) if major_smiles else Counter()

    fg_delta = _delta_summary(product_fg, major_fg)
    pg_delta = _delta_summary(product_pg, major_pg)
    halogen_delta = _delta_summary(product_halogen_atoms, major_halogen_atoms)
    combined_fg_delta = _delta_summary(product_fg, combined_precursor_fg)
    combined_pg_delta = _delta_summary(product_pg, combined_precursor_pg)
    combined_halogen_delta = _delta_summary(product_halogen_atoms, combined_precursor_halogen_atoms)

    delta_codes: List[str] = []
    findings: List[Dict[str, Any]] = []

    if pg_delta["increased"]:
        delta_codes.append("protecting_group_added")
        findings.append(make_finding(
            code="protecting_group_added",
            severity="info",
            source="fg_delta",
            message="protecting-group-like functionality increased in the product",
            evidence={
                "protecting_group_increased": pg_delta["increased"],
                "fg_increased": {},
            },
        ))
    if pg_delta["decreased"]:
        delta_codes.append("protecting_group_removed")
        findings.append(make_finding(
            code="protecting_group_removed",
            severity="info",
            source="fg_delta",
            message="protecting-group-like functionality decreased in the product",
            evidence={
                "protecting_group_decreased": pg_delta["decreased"],
                "fg_decreased": {},
            },
        ))

    if halogen_delta["increased"]:
        delta_codes.append("halogen_added")
        findings.append(make_finding(
            code="halogen_added",
            severity="info",
            source="fg_delta",
            message="halogen atom count increased in the product",
            evidence={"halogen_increased": halogen_delta["increased"]},
        ))
    if halogen_delta["decreased"]:
        delta_codes.append("halogen_removed")
        findings.append(make_finding(
            code="halogen_removed",
            severity="info",
            source="fg_delta",
            message="halogen atom count decreased in the product",
            evidence={"halogen_decreased": halogen_delta["decreased"]},
        ))

    precursor_nitro = _sum_for(major_fg, _NITRO_GROUPS)
    product_nitro = _sum_for(product_fg, _NITRO_GROUPS)
    precursor_amine = _sum_for(major_fg, _AMINE_GROUPS)
    product_amine = _sum_for(product_fg, _AMINE_GROUPS)
    if precursor_nitro > product_nitro and product_amine > precursor_amine:
        delta_codes.append("nitro_to_amine")
        findings.append(make_finding(
            code="nitro_to_amine",
            severity="info",
            source="fg_delta",
            message="nitro functionality decreased while amine functionality increased",
            evidence={
                "precursor_nitro_count": precursor_nitro,
                "product_nitro_count": product_nitro,
                "precursor_amine_count": precursor_amine,
                "product_amine_count": product_amine,
            },
        ))

    precursor_carbonyl = _sum_for(major_fg, _CARBONYL_GROUPS)
    product_carbonyl = _sum_for(product_fg, _CARBONYL_GROUPS)
    precursor_alcohol = _sum_for(major_fg, _ALCOHOL_GROUPS)
    product_alcohol = _sum_for(product_fg, _ALCOHOL_GROUPS)
    if precursor_carbonyl > product_carbonyl and product_alcohol > precursor_alcohol:
        delta_codes.append("carbonyl_to_alcohol")
        findings.append(make_finding(
            code="carbonyl_to_alcohol",
            severity="info",
            source="fg_delta",
            message="carbonyl functionality decreased while alcohol functionality increased",
            evidence={
                "precursor_carbonyl_count": precursor_carbonyl,
                "product_carbonyl_count": product_carbonyl,
                "precursor_alcohol_count": precursor_alcohol,
                "product_alcohol_count": product_alcohol,
            },
        ))
    if product_carbonyl > precursor_carbonyl and precursor_alcohol > product_alcohol:
        delta_codes.append("alcohol_to_carbonyl")
        findings.append(make_finding(
            code="alcohol_to_carbonyl",
            severity="info",
            source="fg_delta",
            message="alcohol functionality decreased while carbonyl functionality increased",
            evidence={
                "precursor_carbonyl_count": precursor_carbonyl,
                "product_carbonyl_count": product_carbonyl,
                "precursor_alcohol_count": precursor_alcohol,
                "product_alcohol_count": product_alcohol,
            },
        ))

    precursor_halide_fg = _sum_for(major_fg, _HALIDE_GROUPS)
    product_halide_fg = _sum_for(product_fg, _HALIDE_GROUPS)
    if product_halide_fg > precursor_halide_fg and precursor_alcohol > product_alcohol:
        delta_codes.append("alcohol_to_halide")
        findings.append(make_finding(
            code="alcohol_to_halide",
            severity="info",
            source="fg_delta",
            message="alcohol functionality decreased while organic halide functionality increased",
            evidence={
                "precursor_alcohol_count": precursor_alcohol,
                "product_alcohol_count": product_alcohol,
                "precursor_halide_fg_count": precursor_halide_fg,
                "product_halide_fg_count": product_halide_fg,
            },
        ))
    if precursor_halide_fg > product_halide_fg and product_alcohol > precursor_alcohol:
        delta_codes.append("halide_to_alcohol")
        findings.append(make_finding(
            code="halide_to_alcohol",
            severity="info",
            source="fg_delta",
            message="organic halide functionality decreased while alcohol functionality increased",
            evidence={
                "precursor_alcohol_count": precursor_alcohol,
                "product_alcohol_count": product_alcohol,
                "precursor_halide_fg_count": precursor_halide_fg,
                "product_halide_fg_count": product_halide_fg,
            },
        ))

    if fg_delta["increased"] or fg_delta["decreased"]:
        delta_codes.append("fg_count_change")

    delta_codes = list(dict.fromkeys(delta_codes))
    action_context = dict(action_context or {})

    return {
        "ok": True,
        "applies": bool(delta_codes),
        "source": "fg_delta",
        "risk_level": "low" if delta_codes else "none",
        "reaction_category": reaction_category or "",
        "product": canonical(product_smiles) or product_smiles,
        "precursors": [canonical(smi) or smi for smi in precursor_smiles],
        "major_precursor_index": major.get("index"),
        "major_precursor": major_smiles,
        "major_precursor_mcs_coverage": major.get("mcs_coverage", 0.0),
        "fg_delta": fg_delta,
        "combined_fg_delta": combined_fg_delta,
        "protecting_group_delta": pg_delta,
        "combined_protecting_group_delta": combined_pg_delta,
        "halogen_atom_delta": halogen_delta,
        "combined_halogen_atom_delta": combined_halogen_delta,
        "delta_codes": delta_codes,
        "findings": findings,
        "violations": [],
        "action_context": action_context,
        "summary": (
            "no functional-group delta"
            if not delta_codes
            else "functional-group delta: " + ", ".join(delta_codes)
        ),
    }
