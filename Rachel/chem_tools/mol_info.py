"""M1: 分子信息与图结构 — 分析分子静态属性（原子、键、环、立体化学、骨架）。"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from ._rdkit_utils import canonical, load_template, parse_mol, smarts_match, tanimoto

# ---------------------------------------------------------------------------
# Bond type name mapping
# ---------------------------------------------------------------------------

_BOND_TYPE_NAME = {
    Chem.BondType.SINGLE: "SINGLE",
    Chem.BondType.DOUBLE: "DOUBLE",
    Chem.BondType.TRIPLE: "TRIPLE",
    Chem.BondType.AROMATIC: "AROMATIC",
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_atom_info(mol: Chem.Mol, idx: int) -> Dict[str, Any]:
    """Return detailed info for a single atom."""
    atom = mol.GetAtomWithIdx(idx)
    return {
        "idx": idx,
        "element": atom.GetSymbol(),
        "charge": atom.GetFormalCharge(),
        "hybridization": str(atom.GetHybridization()),
        "in_ring": atom.IsInRing(),
        "aromatic": atom.GetIsAromatic(),
        "num_hs": atom.GetTotalNumHs(),
        "neighbors": sorted(n.GetIdx() for n in atom.GetNeighbors()),
    }


def _get_bond_info(mol: Chem.Mol, bond: Chem.Bond) -> Dict[str, Any]:
    """Return detailed info for a single bond."""
    return {
        "idx": bond.GetIdx(),
        "atoms": [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()],
        "bond_type": _BOND_TYPE_NAME.get(bond.GetBondType(), str(bond.GetBondType())),
        "in_ring": bond.IsInRing(),
        "aromatic": bond.GetIsAromatic(),
        "conjugated": bond.GetIsConjugated(),
    }


def _analyze_rings(mol: Chem.Mol) -> List[Dict[str, Any]]:
    """Return ring info list from mol.GetRingInfo()."""
    ri = mol.GetRingInfo()
    rings: List[Dict[str, Any]] = []
    for ring_atoms in ri.AtomRings():
        is_aromatic = all(
            mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms
        )
        is_heterocyclic = any(
            mol.GetAtomWithIdx(a).GetSymbol() != "C" for a in ring_atoms
        )
        rings.append({
            "atoms": list(ring_atoms),
            "size": len(ring_atoms),
            "aromatic": is_aromatic,
            "heterocyclic": is_heterocyclic,
        })
    return rings


def _analyze_stereo(mol: Chem.Mol) -> Dict[str, Any]:
    """Analyze stereochemistry: chiral centers and double-bond stereo."""
    # Chiral centers
    chiral_centers: List[Dict[str, Any]] = []
    try:
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        for atom_idx, label in centers:
            chiral_centers.append({"atom_idx": atom_idx, "label": label})
    except Exception:
        pass

    # Double-bond stereo via FindPotentialStereo
    double_bond_stereo: List[Dict[str, Any]] = []
    try:
        stereo_info = Chem.FindPotentialStereo(mol)
        for si in stereo_info:
            if si.type == Chem.StereoType.Bond_Double:
                double_bond_stereo.append({
                    "bond_idx": si.centeredOn,
                    "descriptor": str(si.descriptor),
                    "specified": str(si.specified),
                })
    except Exception:
        pass

    return {
        "chiral_centers": chiral_centers,
        "double_bond_stereo": double_bond_stereo,
    }


def _analyze_symmetry(mol: Chem.Mol) -> Dict[str, Any]:
    """Approximate symmetry analysis via canonical rank equivalence classes."""
    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy == 0:
        return {"ratio": 1.0, "level": "none"}

    # breakTies=False → topologically equivalent atoms share the same rank
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

    # Count unique ranks among heavy atoms only
    rank_groups: Dict[int, List[int]] = defaultdict(list)
    for idx, rank in enumerate(ranks):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 1:  # skip explicit H
            continue
        rank_groups[rank].append(idx)

    n_unique = len(rank_groups)
    ratio = round(n_unique / n_heavy, 4)

    if ratio >= 0.95:
        level = "none"
    elif ratio >= 0.8:
        level = "low"
    elif ratio >= 0.5:
        level = "moderate"
    else:
        level = "high"

    return {"ratio": ratio, "level": level}


def _analyze_scaffold(mol: Chem.Mol) -> Dict[str, Any]:
    """Analyze Murcko scaffold and ring topology tags.

    Tags (non-exclusive): fused_bicyclic, fused_polycyclic, bridged, spiro,
    macrocycle, strained_small_ring, simple_monocyclic, polycyclic_separated,
    acyclic.
    """
    # Murcko scaffold
    murcko = ""
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        murcko = Chem.MolToSmiles(core, canonical=True)
    except Exception:
        pass

    ri = mol.GetRingInfo()
    atom_rings = [set(r) for r in ri.AtomRings()]
    n_rings = len(atom_rings)

    # Pairwise ring relationships
    fused_pairs: List[Tuple[int, int]] = []
    bridged_pairs: List[Tuple[int, int]] = []
    spiro_pairs: List[Tuple[int, int]] = []

    for i in range(n_rings):
        for j in range(i + 1, n_rings):
            shared = atom_rings[i] & atom_rings[j]
            ns = len(shared)
            if ns >= 3:
                bridged_pairs.append((i, j))
            elif ns == 2:
                fused_pairs.append((i, j))
            elif ns == 1:
                spiro_pairs.append((i, j))

    # Fused ring clusters via DFS
    fused_clusters = _cluster_rings(fused_pairs)

    # Tags
    tags: List[str] = []

    # Macrocycle: any ring >= 12 atoms
    ring_sizes = [len(r) for r in atom_rings]
    if any(s >= 12 for s in ring_sizes):
        tags.append("macrocycle")

    # Strained small ring: 3 or 4 membered
    if any(s <= 4 for s in ring_sizes):
        tags.append("strained_small_ring")

    # Bridged
    if bridged_pairs:
        tags.append("bridged")

    # Spiro
    if spiro_pairs:
        tags.append("spiro")

    # Fused
    if fused_clusters:
        max_cluster = max(len(c) for c in fused_clusters)
        if max_cluster >= 3:
            tags.append("fused_polycyclic")
        else:
            tags.append("fused_bicyclic")

    # Separated polycyclic
    if n_rings >= 2 and not tags:
        tags.append("polycyclic_separated")

    # Simple monocyclic
    if n_rings == 1 and not tags:
        tags.append("simple_monocyclic")

    # Acyclic
    if n_rings == 0:
        tags.append("acyclic")

    return {
        "murcko": murcko,
        "tags": tags,
        "fused_pairs": len(fused_pairs),
        "bridged_pairs": len(bridged_pairs),
        "spiro_pairs": len(spiro_pairs),
    }


def _cluster_rings(pairs: List[Tuple[int, int]]) -> List[List[int]]:
    """Cluster ring indices into connected components via DFS."""
    if not pairs:
        return []
    adj: Dict[int, Set[int]] = {}
    for a, b in pairs:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    visited: Set[int] = set()
    clusters: List[List[int]] = []
    for node in adj:
        if node in visited:
            continue
        cluster: List[int] = []
        stack = [node]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            cluster.append(cur)
            for nb in adj.get(cur, set()):
                if nb not in visited:
                    stack.append(nb)
        cluster.sort()
        clusters.append(cluster)
    return clusters


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyze_molecule(smiles: str) -> Dict[str, Any]:
    """Comprehensive molecule analysis — main entry point for M1.

    Returns a dict with keys: smiles, formula, mw, descriptors, atoms, bonds,
    rings, stereo, symmetry, scaffold.

    Invalid SMILES returns ``{"ok": False, "error": "...", "input": smiles}``.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    can_smi = canonical(smiles) or smiles

    # Molecular formula & weight
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = round(Descriptors.MolWt(mol), 2)

    # Descriptors
    descriptors = {
        "logP": round(Descriptors.MolLogP(mol), 2),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heterocycle_count": rdMolDescriptors.CalcNumHeterocycles(mol),
    }

    # Atoms
    atoms = [_get_atom_info(mol, idx) for idx in range(mol.GetNumAtoms())]

    # Bonds
    bonds = [_get_bond_info(mol, bond) for bond in mol.GetBonds()]

    # Rings
    rings = _analyze_rings(mol)

    # Stereo
    stereo = _analyze_stereo(mol)

    # Symmetry
    symmetry = _analyze_symmetry(mol)

    # Scaffold
    scaffold = _analyze_scaffold(mol)

    return {
        "smiles": can_smi,
        "formula": formula,
        "mw": mw,
        "descriptors": descriptors,
        "atoms": atoms,
        "bonds": bonds,
        "rings": rings,
        "stereo": stereo,
        "symmetry": symmetry,
        "scaffold": scaffold,
    }


# ---------------------------------------------------------------------------
# Known scaffold matching
# ---------------------------------------------------------------------------


def _murcko_smiles(mol: Chem.Mol) -> str:
    """Return canonical Murcko scaffold SMILES, or empty string on failure."""
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        if core is None or core.GetNumAtoms() == 0:
            return ""
        return Chem.MolToSmiles(core, canonical=True)
    except Exception:
        return ""


def _score_confidence(method: str, coverage: float, priority: int) -> Dict[str, Any]:
    """Compute heuristic confidence for a scaffold match.

    Returns ``{"value": float, "is_heuristic": True}``.
    Weights are hand-tuned — the ``is_heuristic`` flag makes this explicit.
    """
    if method == "smarts":
        score = 0.55 + 0.35 * coverage + min(priority, 5) * 0.02
    elif method == "murcko_exact":
        score = 0.5 + min(priority, 5) * 0.03
    else:  # murcko_substructure
        score = 0.4 + 0.25 * coverage + min(priority, 5) * 0.02
    return {
        "value": round(min(score, 0.99), 3),
        "is_heuristic": True,
    }


def match_known_scaffolds(
    smiles: str,
    max_hits: int = 5,
) -> Dict[str, Any]:
    """Match a molecule against the known scaffold library.

    Loads ``templates/known_scaffolds.json`` via :func:`load_template`.
    For each scaffold entry the function first attempts a SMARTS substructure
    match (``method="smarts"``).  If that fails it falls back to Murcko
    scaffold comparison (exact match → ``"murcko_exact"``, substructure →
    ``"murcko_substructure"``).

    Returns
    -------
    dict
        ``{"ok": True, "hits": [...]}`` on success.
        ``{"ok": False, "error": "...", "input": smiles}`` for invalid SMILES.

    Each hit contains: scaffold_id, name, class, method, coverage, confidence,
    matched_atoms, reference_route_tags.  ``confidence`` is annotated with
    ``is_heuristic``.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    n_heavy = mol.GetNumHeavyAtoms()
    target_murcko = _murcko_smiles(mol)

    # Load scaffold library
    raw_data = load_template("known_scaffolds.json")
    scaffolds: List[Dict[str, Any]]
    if isinstance(raw_data, dict):
        scaffolds = raw_data.get("scaffolds", [])
    elif isinstance(raw_data, list):
        scaffolds = raw_data
    else:
        scaffolds = []

    hits: List[Dict[str, Any]] = []

    for entry in scaffolds:
        smarts_str = entry.get("smarts", "")
        murcko_smi = entry.get("murcko_smiles", "")
        priority = int(entry.get("priority", 1))

        method = ""
        matched_atoms: Tuple[int, ...] = ()
        coverage = 0.0

        # --- 1. Try SMARTS matching first ---
        if smarts_str:
            matches = smarts_match(mol, smarts_str)
            if matches:
                # Pick the match with the most unique atoms
                best = max(matches, key=lambda m: len(set(m)))
                matched_atoms = best
                coverage = len(set(best)) / max(n_heavy, 1)
                method = "smarts"

        # --- 2. Murcko fallback ---
        if not method and murcko_smi and target_murcko:
            # Exact Murcko match
            if target_murcko == Chem.MolToSmiles(
                Chem.MolFromSmiles(murcko_smi), canonical=True
            ) if Chem.MolFromSmiles(murcko_smi) else False:
                method = "murcko_exact"
                coverage = 1.0
                murcko_mol = parse_mol(murcko_smi)
                if murcko_mol is not None:
                    sub_matches = smarts_match(mol, Chem.MolToSmarts(murcko_mol))
                    if sub_matches:
                        matched_atoms = max(sub_matches, key=lambda m: len(set(m)))
            else:
                # Murcko substructure match
                murcko_mol = parse_mol(murcko_smi)
                if murcko_mol is not None:
                    sub_matches = smarts_match(mol, Chem.MolToSmarts(murcko_mol))
                    if sub_matches:
                        best = max(sub_matches, key=lambda m: len(set(m)))
                        matched_atoms = best
                        coverage = len(set(best)) / max(n_heavy, 1)
                        method = "murcko_substructure"

        if not method:
            continue

        conf = _score_confidence(method, coverage, priority)
        hits.append({
            "scaffold_id": entry.get("id", ""),
            "name": entry.get("name", ""),
            "class": entry.get("class", "unknown"),
            "method": method,
            "coverage": round(coverage, 3),
            "confidence": conf["value"],
            "is_heuristic": conf["is_heuristic"],
            "matched_atoms": sorted(set(int(i) for i in matched_atoms)),
            "reference_route_tags": list(entry.get("reference_route_tags", [])),
        })

    # Sort by confidence descending
    hits.sort(key=lambda h: -h["confidence"])
    top_hits = hits[:max(1, max_hits)]

    return {
        "ok": True,
        "hits": top_hits,
    }
