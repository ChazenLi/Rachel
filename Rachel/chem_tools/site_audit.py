from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

from ._rdkit_utils import canonical, parse_mol
from .mol_info import analyze_molecule


_SITE_RETENTION_REASONING_KEYWORDS = (
    "same site",
    "same carbon",
    "same atom",
    "same c2",
    "same c3",
    "site-retentive",
    "site retentive",
    "same-site",
    "同位点",
    "同一位点",
    "同一碳位点",
    "同一碳",
    "只改",
    "保持不变",
)


def reasoning_claims_site_retention(reasoning: str) -> bool:
    """Return True when reasoning explicitly claims same-site rollback/retention."""
    text = (reasoning or "").strip().lower()
    if not text:
        return False
    return any(keyword in text for keyword in _SITE_RETENTION_REASONING_KEYWORDS)


def _murcko_smiles(mol: Optional[Chem.Mol]) -> str:
    if mol is None:
        return ""
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        if core is None or core.GetNumAtoms() == 0:
            return ""
        return Chem.MolToSmiles(core, canonical=True)
    except Exception:
        return ""


def _largest_ring_core(mol: Optional[Chem.Mol]) -> Tuple[str, Optional[Chem.Mol]]:
    """Return the canonical fragment for the largest connected ring system."""
    if mol is None:
        return "", None

    ring_info = mol.GetRingInfo()
    atom_rings = [set(ring) for ring in ring_info.AtomRings()]
    if not atom_rings:
        return "", None

    clusters: List[set[int]] = []
    for ring_atoms in atom_rings:
        attached_clusters = [cluster for cluster in clusters if cluster & ring_atoms]
        if not attached_clusters:
            clusters.append(set(ring_atoms))
            continue

        merged = set(ring_atoms)
        for cluster in attached_clusters:
            merged.update(cluster)
            clusters.remove(cluster)
        clusters.append(merged)

    def _cluster_score(atom_ids: set[int]) -> Tuple[int, int, int]:
        hetero = sum(
            1 for atom_idx in atom_ids
            if mol.GetAtomWithIdx(atom_idx).GetSymbol() != "C"
        )
        aromatic = sum(
            1 for atom_idx in atom_ids
            if mol.GetAtomWithIdx(atom_idx).GetIsAromatic()
        )
        return (len(atom_ids), hetero, aromatic)

    largest_cluster = max(clusters, key=_cluster_score)
    core_atom_ids = sorted(largest_cluster)
    core_smiles = Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=core_atom_ids,
        canonical=True,
    )
    try:
        core_smarts = Chem.MolFragmentToSmarts(mol, atomsToUse=core_atom_ids)
    except Exception:
        core_smarts = ""
    core_query = Chem.MolFromSmarts(core_smarts) if core_smarts else None
    return core_smiles, core_query


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


def _normalize_attachment_smiles(smiles: str) -> str:
    """Collapse RDKit dummy labels so audit output stays stable and compact."""
    return re.sub(r"\[\d+\*\]", "[*]", smiles or "")


def _attachment_fragment_smiles(
    mol: Chem.Mol,
    bond_idx: int,
    external_atom_idx: int,
) -> str:
    try:
        frag_mol = Chem.FragmentOnBonds(Chem.Mol(mol), [bond_idx], addDummies=True)
        for atom_ids in Chem.GetMolFrags(frag_mol):
            if external_atom_idx in atom_ids:
                smiles = Chem.MolFragmentToSmiles(
                    frag_mol,
                    atomsToUse=sorted(atom_ids),
                    canonical=True,
                )
                return _normalize_attachment_smiles(smiles)
    except Exception:
        pass
    return ""


def _collect_site_substituents(
    mol: Chem.Mol,
    scaffold_mol: Chem.Mol,
    scaffold_match: Sequence[int],
) -> Dict[int, Dict[str, Any]]:
    scaffold_atom_set = set(int(idx) for idx in scaffold_match)
    site_map: Dict[int, Dict[str, Any]] = {}

    for scaffold_idx, mol_atom_idx in enumerate(scaffold_match):
        mol_atom = mol.GetAtomWithIdx(int(mol_atom_idx))
        scaffold_atom = scaffold_mol.GetAtomWithIdx(int(scaffold_idx))
        descriptors: List[str] = []
        for neighbor in mol_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in scaffold_atom_set:
                continue
            bond = mol.GetBondBetweenAtoms(int(mol_atom_idx), int(neighbor_idx))
            if bond is None:
                continue
            descriptor = _attachment_fragment_smiles(mol, bond.GetIdx(), neighbor_idx)
            if descriptor:
                descriptors.append(descriptor)

        if descriptors:
            site_map[int(scaffold_idx)] = {
                "scaffold_idx": int(scaffold_idx),
                "symbol": scaffold_atom.GetSymbol(),
                "mol_atom_idx": int(mol_atom_idx),
                "substituents": sorted(descriptors),
            }

    return site_map


def _compare_site_maps(
    scaffold_mol: Chem.Mol,
    product_map: Dict[int, Dict[str, Any]],
    precursor_map: Dict[int, Dict[str, Any]],
) -> Dict[str, Any]:
    changed_sites: List[Dict[str, Any]] = []
    preserved_sites: List[Dict[str, Any]] = []
    site_rows: List[Dict[str, Any]] = []

    for scaffold_idx in sorted(set(product_map) | set(precursor_map)):
        product_entry = product_map.get(scaffold_idx, {})
        precursor_entry = precursor_map.get(scaffold_idx, {})
        product_subs = list(product_entry.get("substituents", []))
        precursor_subs = list(precursor_entry.get("substituents", []))
        scaffold_atom = scaffold_mol.GetAtomWithIdx(int(scaffold_idx))
        site_label = (
            f"scaf_{scaffold_idx}:{scaffold_atom.GetSymbol()}"
            f"(product_atom={product_entry.get('mol_atom_idx', -1)},"
            f" precursor_atom={precursor_entry.get('mol_atom_idx', -1)})"
        )
        site_row = {
            "site_label": site_label,
            "scaffold_idx": int(scaffold_idx),
            "product_substituents": product_subs,
            "precursor_substituents": precursor_subs,
        }
        if product_subs == precursor_subs:
            site_row["status"] = "preserved"
            preserved_sites.append(site_row)
        else:
            site_row["status"] = "changed"
            changed_sites.append(site_row)
        site_rows.append(site_row)

    return {
        "site_rows": site_rows,
        "changed_sites": changed_sites,
        "preserved_sites": preserved_sites,
    }


def _summarize_site_audit(
    *,
    site_retentive: bool,
    pass_audit: bool,
    changed_sites: Sequence[Dict[str, Any]],
    preserved_sites: Sequence[Dict[str, Any]],
) -> str:
    if not site_retentive:
        return "No same-core major precursor was identified, so site-retention audit did not apply."

    if pass_audit:
        if changed_sites:
            changed = changed_sites[0]
            return (
                "Site-retentive audit passed: only "
                f"{changed['site_label']} changed from {changed['precursor_substituents']} "
                f"to {changed['product_substituents']}; preserved sites="
                f"{[row['site_label'] for row in preserved_sites]}"
            )
        return (
            "Site-retentive audit passed: all detected scaffold attachment sites were preserved; "
            f"preserved sites={[row['site_label'] for row in preserved_sites]}"
        )

    return (
        "Site anchor drift detected: multiple scaffold attachment sites changed under a same-core "
        f"proposal; changed sites={[row['site_label'] for row in changed_sites]}"
    )


def audit_site_retention(product_smiles: str, precursors: Sequence[str]) -> Dict[str, Any]:
    """Audit whether a proposed major precursor preserves the same scaffold attachment sites.

    This is intentionally read-only and compact. It keeps the full atom/bond comparison
    in the chemistry tool layer and exports only the minimal gate-friendly summary that
    the sandbox/session flow needs.
    """
    product_can = canonical(product_smiles) or product_smiles
    product_mol = parse_mol(product_can)
    if product_mol is None:
        return {"ok": False, "error": "invalid product smiles", "input": product_smiles}

    product_info = analyze_molecule(product_can)
    product_murcko = product_info.get("scaffold", {}).get("murcko", "")
    product_ring_core, product_ring_core_mol = _largest_ring_core(product_mol)
    product_tags = list(product_info.get("scaffold", {}).get("tags", []))

    candidates: List[Dict[str, Any]] = []
    for precursor_idx, precursor_smiles in enumerate(precursors):
        precursor_can = canonical(precursor_smiles) or precursor_smiles
        precursor_mol = parse_mol(precursor_can)
        if precursor_mol is None:
            continue

        precursor_info = analyze_molecule(precursor_can)
        precursor_murcko = precursor_info.get("scaffold", {}).get("murcko", "")
        precursor_ring_core, _ = _largest_ring_core(precursor_mol)
        same_ring_core = bool(product_ring_core and product_ring_core == precursor_ring_core)
        murcko_equal = bool(product_murcko and product_murcko == precursor_murcko)
        coverage = _mcs_coverage(product_mol, precursor_mol)
        score = (
            (2000.0 if same_ring_core else 0.0)
            + (1000.0 if murcko_equal else 0.0)
            + coverage * 100.0
            + precursor_mol.GetNumHeavyAtoms()
        )
        candidates.append({
            "precursor_index": precursor_idx,
            "precursor": precursor_can,
            "precursor_mol": precursor_mol,
            "precursor_info": precursor_info,
            "same_ring_core": same_ring_core,
            "murcko_equal": murcko_equal,
            "mcs_coverage": coverage,
            "score": score,
        })

    if not candidates:
        return {
            "ok": True,
            "site_retentive": False,
            "strict_audit_required": False,
            "pass": True,
            "summary": "No valid precursor candidate was available for same-core site audit.",
        }

    candidates.sort(key=lambda item: (-item["score"], item["precursor_index"]))
    best = candidates[0]
    major_precursor_mol = best["precursor_mol"]
    major_precursor = best["precursor"]
    same_ring_core = bool(best.get("same_ring_core"))
    same_murcko = bool(best["murcko_equal"])

    # Historical candidate gate kept here as a comment instead of deleting it.
    # Reason: the first draft treated any high-MCS precursor as site-retentive,
    # which was too permissive for the user's fused heteroaryl error class.
    # site_retentive = bool(best["mcs_coverage"] >= 0.75)
    site_retentive = same_ring_core or same_murcko

    strict_audit_required = bool(
        site_retentive and (
            any(tag in product_tags for tag in ("fused_bicyclic", "fused_polycyclic", "bridged", "spiro"))
            or (
                product_info.get("descriptors", {}).get("aromatic_ring_count", 0) >= 1
                and product_mol.GetNumHeavyAtoms() >= 10
            )
        )
    )

    if not site_retentive:
        return {
            "ok": True,
            "site_retentive": False,
            "strict_audit_required": strict_audit_required,
            "major_precursor_index": best["precursor_index"],
            "major_precursor": major_precursor,
            "same_ring_core": same_ring_core,
            "murcko_equal": same_murcko,
            "mcs_coverage": best["mcs_coverage"],
            "pass": True,
            "summary": "No same-core Murcko-preserving major precursor was detected, so site-retention gate does not apply.",
        }

    # Historical core choice kept as a comment instead of deleting it.
    # Reason: pure Murcko was too strict for same-indole-core / wrong-site swaps,
    # because moving a large appended aryl-ethynyl group can change the Murcko string.
    # scaffold_mol = MurckoScaffold.GetScaffoldForMol(product_mol)
    scaffold_mol = product_ring_core_mol if same_ring_core and product_ring_core_mol is not None else MurckoScaffold.GetScaffoldForMol(product_mol)
    if scaffold_mol is None or scaffold_mol.GetNumAtoms() == 0:
        return {
            "ok": True,
            "site_retentive": True,
            "strict_audit_required": strict_audit_required,
            "major_precursor_index": best["precursor_index"],
            "major_precursor": major_precursor,
            "same_ring_core": same_ring_core,
            "murcko_equal": same_murcko,
            "mcs_coverage": best["mcs_coverage"],
            "pass": False,
            "summary": "Same-core precursor detected but Murcko scaffold expansion failed, so explicit site audit is missing.",
        }

    product_matches = list(product_mol.GetSubstructMatches(scaffold_mol, uniquify=False))[:6]
    precursor_matches = list(major_precursor_mol.GetSubstructMatches(scaffold_mol, uniquify=False))[:6]
    if not product_matches or not precursor_matches:
        return {
            "ok": True,
            "site_retentive": True,
            "strict_audit_required": strict_audit_required,
            "major_precursor_index": best["precursor_index"],
            "major_precursor": major_precursor,
            "same_ring_core": same_ring_core,
            "murcko_equal": same_murcko,
            "mcs_coverage": best["mcs_coverage"],
            "pass": False,
            "summary": "Same-core precursor detected but scaffold atom mapping failed, so explicit site audit is missing.",
        }

    best_comparison: Optional[Dict[str, Any]] = None
    for product_match in product_matches:
        product_site_map = _collect_site_substituents(product_mol, scaffold_mol, product_match)
        for precursor_match in precursor_matches:
            precursor_site_map = _collect_site_substituents(major_precursor_mol, scaffold_mol, precursor_match)
            comparison = _compare_site_maps(scaffold_mol, product_site_map, precursor_site_map)
            comparison["product_site_map"] = product_site_map
            comparison["precursor_site_map"] = precursor_site_map
            comparison["product_match"] = list(product_match)
            comparison["precursor_match"] = list(precursor_match)
            comparison["changed_site_count"] = len(comparison["changed_sites"])
            comparison["preserved_site_count"] = len(comparison["preserved_sites"])

            if best_comparison is None:
                best_comparison = comparison
                continue

            if comparison["changed_site_count"] < best_comparison["changed_site_count"]:
                best_comparison = comparison
                continue

            if (
                comparison["changed_site_count"] == best_comparison["changed_site_count"]
                and comparison["preserved_site_count"] > best_comparison["preserved_site_count"]
            ):
                best_comparison = comparison

    if best_comparison is None:
        return {
            "ok": True,
            "site_retentive": True,
            "strict_audit_required": strict_audit_required,
            "major_precursor_index": best["precursor_index"],
            "major_precursor": major_precursor,
            "same_ring_core": same_ring_core,
            "murcko_equal": same_murcko,
            "mcs_coverage": best["mcs_coverage"],
            "pass": False,
            "summary": "Same-core precursor detected but no scaffold-site comparison could be completed.",
        }

    changed_sites = list(best_comparison["changed_sites"])
    preserved_sites = list(best_comparison["preserved_sites"])
    changed_site_count = int(best_comparison["changed_site_count"])
    pass_audit = changed_site_count <= 1

    return {
        "ok": True,
        "site_retentive": True,
        "strict_audit_required": strict_audit_required,
        "major_precursor_index": best["precursor_index"],
        "major_precursor": major_precursor,
        "same_ring_core": same_ring_core,
        "murcko_equal": same_murcko,
        "mcs_coverage": best["mcs_coverage"],
        "changed_site_count": changed_site_count,
        "changed_sites": changed_sites,
        "preserved_sites": preserved_sites,
        "site_rows": list(best_comparison["site_rows"]),
        "pass": pass_audit,
        "summary": _summarize_site_audit(
            site_retentive=True,
            pass_audit=pass_audit,
            changed_sites=changed_sites,
            preserved_sites=preserved_sites,
        ),
    }
