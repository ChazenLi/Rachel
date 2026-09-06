"""M6: CS 合成难度评分 — 自定义Complicated Score评分系统，评估分子合成难度。"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Set

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

from Rachel.knowledge import get_base_profile

from ._rdkit_utils import parse_mol, smarts_match, tanimoto

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CS_DIMENSIONS: Dict[str, float] = {
    "size": 0.55,       # 分子大小（MW）
    "ring": 0.65,       # 环系统拓扑（稠合/桥环/螺环/大环）
    "stereo": 0.55,     # 立体化学负担
    "hetero": 0.40,     # 杂原子密度与多样性
    "symmetry": -0.20,  # 对称性折扣
    "fg_density": 0.35, # 官能团密度
}

CS_TRIVIAL: float = 2.25  # ≤2.25: trivial, is_terminal=true
CS_MODERATE: float = 6.0  # ≤6.0: moderate; >6.0: complex

def _cs_rules(knowledge_profile=None):
    profile = knowledge_profile or get_base_profile()
    payload = profile.get("chem.cs_smarts")
    entries = payload.get("entries", []) if isinstance(payload, dict) else []
    return profile, [entry for entry in entries if isinstance(entry, dict)]


def _knowledge_refs(profile, entry_ids: Set[str]) -> List[Dict[str, str]]:
    return [profile.source("chem.cs_smarts", entry_id) for entry_id in sorted(entry_ids)]


# ---------------------------------------------------------------------------
# Dimension helpers
# ---------------------------------------------------------------------------

def _dim_size(mol: Chem.Mol) -> float:
    """Size dimension based on heavy atom count. Cap 3.0.

    Heavy atom count is a better proxy for synthetic complexity than MW
    because heteroatom-heavy small molecules (e.g. CF3 groups) inflate MW
    without proportionally increasing step count.
    """
    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy <= 6:
        return 0.0
    # Smooth log curve: 10 heavy → 0.35, 20 → 1.05, 30 → 1.50, 50 → 2.10
    return min(math.log2(n_heavy / 6) * 1.2, 3.0)


def _dim_ring(mol: Chem.Mol) -> float:
    """Ring dimension: topology-based scoring. Cap 5.0.

    Each ring adds synthetic steps. Fused/bridged/spiro systems are
    disproportionately harder due to stereochemical and strain constraints.
    """
    ri = mol.GetRingInfo()
    atom_rings = [set(r) for r in ri.AtomRings()]
    n_rings = len(atom_rings)
    if n_rings == 0:
        return 0.0

    n_fused = 0
    n_bridged = 0
    n_spiro = 0
    for i in range(n_rings):
        for j in range(i + 1, n_rings):
            shared = len(atom_rings[i] & atom_rings[j])
            if shared >= 3:
                n_bridged += 1
            elif shared == 2:
                n_fused += 1
            elif shared == 1:
                n_spiro += 1

    n_macrocycle = sum(1 for r in atom_rings if len(r) > 8)
    n_heterocyclic = sum(
        1 for r in atom_rings
        if any(mol.GetAtomWithIdx(idx).GetSymbol() != "C" for idx in r)
    )

    # Independent rings (not part of any fused/bridged/spiro pair)
    n_paired = n_fused + n_bridged + n_spiro
    n_independent = max(n_rings - n_paired, 0)

    score = 0.0
    score += min(n_independent * 0.4, 1.6)   # each independent ring
    score += n_fused * 0.8                     # fused pairs
    score += n_bridged * 1.5                   # bridged pairs (norbornane etc.)
    score += n_spiro * 1.0                     # spiro junctions
    score += n_macrocycle * 1.5                # macrocycles (>8 atoms)
    score += min(n_heterocyclic * 0.2, 1.0)   # heterocyclic rings
    return min(score, 5.0)


def _dim_stereo(mol: Chem.Mol) -> float:
    """Stereo dimension: chiral burden. Cap 3.5."""
    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy == 0:
        return 0.0

    try:
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    except Exception:
        return 0.0

    n_chiral = len(chiral_centers)
    if n_chiral == 0:
        return 0.0

    # Base contribution (diminishing returns)
    score = 0.4 * math.sqrt(n_chiral)

    # Density bonus
    density = n_chiral / n_heavy
    if density > 0.15:
        score += (density - 0.15) * 5.0

    # Consecutive chiral centers (adjacent chiral atoms)
    chiral_idxs = {idx for idx, _ in chiral_centers}
    n_consecutive = 0
    for idx in chiral_idxs:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in chiral_idxs and nbr.GetIdx() > idx:
                n_consecutive += 1
    score += min(n_consecutive * 0.3, 1.5)

    return min(score, 3.5)


def _dim_hetero(mol: Chem.Mol) -> float:
    """Hetero dimension: heteroatom density & variety. Cap 2.5."""
    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy == 0:
        return 0.0

    hetero_types: set = set()
    n_hetero = 0
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in ("C", "H"):
            n_hetero += 1
            hetero_types.add(sym)

    score = 0.0
    density = n_hetero / n_heavy
    if density > 0.1:
        score += (density - 0.1) * 3.0

    score += min(len(hetero_types) * 0.15, 0.8)
    return min(score, 2.5)


def _dim_symmetry(mol: Chem.Mol) -> float:
    """Symmetry dimension: discount for symmetric molecules.

    Returns a *positive* value representing the discount magnitude.
    The caller applies the negative weight from CS_DIMENSIONS.
    Formula: (0.8 - ratio) * 2.0 if ratio < 0.8, else 0.
    """
    n_heavy = mol.GetNumHeavyAtoms()
    if n_heavy == 0:
        return 0.0

    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    # Only count heavy-atom ranks
    heavy_ranks: set = set()
    for idx, rank in enumerate(ranks):
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 1:
            heavy_ranks.add(rank)

    ratio = len(heavy_ranks) / n_heavy
    if ratio < 0.8:
        return (0.8 - ratio) * 2.0
    return 0.0


def _dim_fg_density(
    mol: Chem.Mol,
    pg_rules: List[Dict[str, Any]],
    complexity_rules: List[Dict[str, Any]],
) -> tuple[float, Set[str]]:
    """FG density: count complexity-adding functional groups. Cap 2.5.

    Counts both protecting groups and synthetic-step-adding FGs
    (amides, sulfonamides, carbamates, ureas, etc.) that each typically
    require a dedicated coupling/formation step.
    """
    # Protecting groups (higher weight — they imply protect+deprotect steps)
    n_pg = 0
    matched_ids: Set[str] = set()
    for rule in pg_rules:
        matches = smarts_match(mol, str(rule.get("smarts", "")))
        n_pg += len(matches)
        if matches:
            matched_ids.add(str(rule["id"]))

    complexity_score = 0.0
    for rule in complexity_rules:
        matches = smarts_match(mol, str(rule.get("smarts", "")))
        complexity_score += len(matches) * float(rule.get("weight", 0.25))
        if matches:
            matched_ids.add(str(rule["id"]))

    score = n_pg * 0.35 + complexity_score
    return min(score, 2.5), matched_ids


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_cs_score(smiles: str, knowledge_profile=None) -> Dict[str, Any]:
    """Compute the CS (Complicated Score) for a molecule.

    Evaluates 6 dimensions (size, ring, stereo, hetero, symmetry, fg_density),
    applies weighted sum, and clamps to [1.0, 10.0].

    Returns
    -------
    dict
        On success::

            {
                "cs_score": float,           # [1.0, 10.0]
                "classification": str,       # "trivial" | "moderate" | "complex"
                "is_terminal": bool,         # True if trivial
                "is_heuristic": True,
                "breakdown": {
                    "size": float,
                    "ring": float,
                    "stereo": float,
                    "hetero": float,
                    "symmetry": float,       # negative value (discount)
                    "fg_density": float
                }
            }

        On invalid SMILES::

            {"ok": False, "error": "invalid SMILES", "input": smiles}
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    profile, rules = _cs_rules(knowledge_profile)
    pg_rules = [rule for rule in rules if rule.get("kind") == "protecting_group"]
    complexity_rules = [
        rule for rule in rules if rule.get("kind") == "complexity_group"
    ]

    # Compute raw dimension values
    raw_size = _dim_size(mol)
    raw_ring = _dim_ring(mol)
    raw_stereo = _dim_stereo(mol)
    raw_hetero = _dim_hetero(mol)
    raw_symmetry = _dim_symmetry(mol)   # positive magnitude
    raw_fg, matched_ids = _dim_fg_density(mol, pg_rules, complexity_rules)

    # Weighted contributions
    w = CS_DIMENSIONS
    c_size = w["size"] * raw_size
    c_ring = w["ring"] * raw_ring
    c_stereo = w["stereo"] * raw_stereo
    c_hetero = w["hetero"] * raw_hetero
    c_symmetry = w["symmetry"] * raw_symmetry   # negative weight → negative contribution
    c_fg = w["fg_density"] * raw_fg

    # Base score 1.0 + weighted sum, clamped to [1.0, 10.0]
    raw_score = 1.0 + c_size + c_ring + c_stereo + c_hetero + c_symmetry + c_fg
    cs_score = max(1.0, min(round(raw_score, 4), 10.0))

    # Classification
    if cs_score <= CS_TRIVIAL:
        classification = "trivial"
        is_terminal = True
    elif cs_score <= CS_MODERATE:
        classification = "moderate"
        is_terminal = False
    else:
        classification = "complex"
        is_terminal = False

    return {
        "cs_score": cs_score,
        "classification": classification,
        "is_terminal": is_terminal,
        "is_heuristic": True,
        "knowledge_profile_hash": profile.digest,
        "knowledge_refs": _knowledge_refs(profile, matched_ids),
        "breakdown": {
            "size": round(c_size, 4),
            "ring": round(c_ring, 4),
            "stereo": round(c_stereo, 4),
            "hetero": round(c_hetero, 4),
            "symmetry": round(c_symmetry, 4),
            "fg_density": round(c_fg, 4),
        },
    }


def classify_complexity(smiles: str, knowledge_profile=None) -> Dict[str, Any]:
    """Quick complexity classification without full breakdown.

    Returns
    -------
    dict
        On success::

            {
                "classification": str,   # "trivial" | "moderate" | "complex"
                "is_terminal": bool,
                "is_heuristic": True
            }

        On invalid SMILES::

            {"ok": False, "error": "invalid SMILES", "input": smiles}
    """
    result = compute_cs_score(smiles, knowledge_profile=knowledge_profile)

    # Propagate error
    if "ok" in result and result["ok"] is False:
        return result

    return {
        "classification": result["classification"],
        "is_terminal": result["is_terminal"],
        "is_heuristic": True,
        "knowledge_profile_hash": result["knowledge_profile_hash"],
        "knowledge_refs": result["knowledge_refs"],
    }


# ---------------------------------------------------------------------------
# score_progress helpers
# ---------------------------------------------------------------------------

def _count_pg_matches(
    mol: Chem.Mol, pg_rules: List[Dict[str, Any]]
) -> tuple[int, Set[str]]:
    """Count protecting-group matches and return the matching rule IDs."""
    n = 0
    matched_ids: Set[str] = set()
    for rule in pg_rules:
        matches = smarts_match(mol, str(rule.get("smarts", "")))
        n += len(matches)
        if matches:
            matched_ids.add(str(rule["id"]))
    return n, matched_ids


def _has_substructure_relation(mol_a: Chem.Mol, mol_b: Chem.Mol) -> bool:
    """Return True if *mol_a* is a substructure of *mol_b* or vice-versa."""
    try:
        if mol_b.HasSubstructMatch(mol_a):
            return True
        if mol_a.HasSubstructMatch(mol_b):
            return True
    except Exception:
        pass
    return False


def _fg_conversion_bonus(
    mol_inter: Chem.Mol,
    mol_target: Chem.Mol,
    conversion_rules: List[Dict[str, Any]],
) -> tuple[float, Set[str]]:
    """Return +0.02 if any FG conversion pair matches (intermediate has source,
    target has dest), else 0.0."""
    matches: List[tuple[float, str]] = []
    for rule in conversion_rules:
        if smarts_match(mol_inter, str(rule.get("source_smarts", ""))) and smarts_match(
            mol_target, str(rule.get("target_smarts", ""))
        ):
            matches.append((float(rule.get("bonus", 0.02)), str(rule["id"])))
    if not matches:
        return 0.0, set()
    return max(bonus for bonus, _entry_id in matches), {
        entry_id for _bonus, entry_id in matches
    }


# ---------------------------------------------------------------------------
# Public API — progress scoring
# ---------------------------------------------------------------------------

def score_progress(
    intermediate: str, target: str, knowledge_profile=None
) -> Dict[str, Any]:
    """Evaluate synthesis progress from *intermediate* toward *target*.

    Scoring components:

    * **Tanimoto** (base): Morgan fingerprint similarity
    * **Substructure bonus** (+0.05): if one molecule is a substructure of the other
    * **Deprotection bonus** (+0.03 per PG): if intermediate has protecting groups
      that target doesn't (virtual deprotection would bring it closer)
    * **FG conversion bonus** (+0.02): if intermediate has functional groups that
      could be converted to target's functional groups

    Returns
    -------
    dict
        On success::

            {
                "progress_score": float,       # [0, 1]
                "tanimoto": float,
                "substructure_bonus": float,
                "deprotection_bonus": float,
                "fg_conversion_bonus": float,
                "is_heuristic": True
            }

        On invalid SMILES::

            {"ok": False, "error": "invalid SMILES", "input": <smiles>}
    """
    mol_inter = parse_mol(intermediate)
    if mol_inter is None:
        return {"ok": False, "error": "invalid SMILES", "input": intermediate}

    mol_target = parse_mol(target)
    if mol_target is None:
        return {"ok": False, "error": "invalid SMILES", "input": target}

    profile, rules = _cs_rules(knowledge_profile)
    pg_rules = [rule for rule in rules if rule.get("kind") == "protecting_group"]
    conversion_rules = [rule for rule in rules if rule.get("kind") == "fg_conversion"]

    # 1. Tanimoto base
    tan_score = tanimoto(mol_inter, mol_target)

    # 2. Substructure bonus
    sub_bonus = 0.05 if _has_substructure_relation(mol_inter, mol_target) else 0.0

    # 3. Deprotection bonus: +0.03 per PG that intermediate has but target doesn't
    pg_inter, inter_pg_ids = _count_pg_matches(mol_inter, pg_rules)
    pg_target, target_pg_ids = _count_pg_matches(mol_target, pg_rules)
    extra_pg = max(pg_inter - pg_target, 0)
    deprot_bonus = extra_pg * 0.03

    # 4. FG conversion bonus
    fg_bonus, conversion_ids = _fg_conversion_bonus(
        mol_inter, mol_target, conversion_rules
    )

    progress = min(1.0, tan_score + sub_bonus + deprot_bonus + fg_bonus)

    return {
        "progress_score": round(progress, 6),
        "tanimoto": round(tan_score, 6),
        "substructure_bonus": round(sub_bonus, 6),
        "deprotection_bonus": round(deprot_bonus, 6),
        "fg_conversion_bonus": round(fg_bonus, 6),
        "is_heuristic": True,
        "knowledge_profile_hash": profile.digest,
        "knowledge_refs": _knowledge_refs(
            profile, inter_pg_ids | target_pg_ids | conversion_ids
        ),
    }


def batch_score_progress(
    intermediates: List[str], target: str, knowledge_profile=None
) -> List[Dict[str, Any]]:
    """Batch-evaluate synthesis progress for multiple intermediates.

    Pre-computes target-side data (Mol, fingerprint, PG count, substructure
    pattern) once, then iterates over *intermediates* to avoid redundant work.

    The i-th result is identical to ``score_progress(intermediates[i], target)``.

    Returns
    -------
    list[dict]
        A list of result dicts, one per intermediate.  Each dict has the same
        schema as :func:`score_progress`.
    """
    mol_target = parse_mol(target)
    if mol_target is None:
        return [
            {"ok": False, "error": "invalid SMILES", "input": target}
            for _ in intermediates
        ]

    profile, rules = _cs_rules(knowledge_profile)
    pg_rules = [rule for rule in rules if rule.get("kind") == "protecting_group"]
    conversion_rules = [rule for rule in rules if rule.get("kind") == "fg_conversion"]

    # Pre-compute target-side data
    fp_target = AllChem.GetMorganFingerprintAsBitVect(mol_target, radius=2, nBits=2048)
    pg_target, target_pg_ids = _count_pg_matches(mol_target, pg_rules)

    # Pre-compile FG conversion destination patterns for target
    _fg_dst_hits: List[bool] = []
    for rule in conversion_rules:
        _fg_dst_hits.append(
            bool(smarts_match(mol_target, str(rule.get("target_smarts", ""))))
        )

    results: List[Dict[str, Any]] = []
    for smi in intermediates:
        mol_inter = parse_mol(smi)
        if mol_inter is None:
            results.append({"ok": False, "error": "invalid SMILES", "input": smi})
            continue

        # 1. Tanimoto (use pre-computed target fp)
        fp_inter = AllChem.GetMorganFingerprintAsBitVect(mol_inter, radius=2, nBits=2048)
        tan_score = DataStructs.TanimotoSimilarity(fp_inter, fp_target)

        # 2. Substructure bonus
        sub_bonus = 0.05 if _has_substructure_relation(mol_inter, mol_target) else 0.0

        # 3. Deprotection bonus
        pg_inter, inter_pg_ids = _count_pg_matches(mol_inter, pg_rules)
        extra_pg = max(pg_inter - pg_target, 0)
        deprot_bonus = extra_pg * 0.03

        # 4. FG conversion bonus (use pre-computed target dst hits)
        conversion_matches: List[tuple[float, str]] = []
        for idx, rule in enumerate(conversion_rules):
            if _fg_dst_hits[idx] and smarts_match(
                mol_inter, str(rule.get("source_smarts", ""))
            ):
                conversion_matches.append(
                    (float(rule.get("bonus", 0.02)), str(rule["id"]))
                )
        fg_bonus = max(
            (bonus for bonus, _entry_id in conversion_matches), default=0.0
        )
        conversion_ids = {entry_id for _bonus, entry_id in conversion_matches}

        progress = min(1.0, tan_score + sub_bonus + deprot_bonus + fg_bonus)

        results.append({
            "progress_score": round(progress, 6),
            "tanimoto": round(tan_score, 6),
            "substructure_bonus": round(sub_bonus, 6),
            "deprotection_bonus": round(deprot_bonus, 6),
            "fg_conversion_bonus": round(fg_bonus, 6),
            "is_heuristic": True,
            "knowledge_profile_hash": profile.digest,
            "knowledge_refs": _knowledge_refs(
                profile, inter_pg_ids | target_pg_ids | conversion_ids
            ),
        })

    return results
