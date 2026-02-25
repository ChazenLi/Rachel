"""M6: CS 合成难度评分 — 自定义Complicated Score评分系统，评估分子合成难度。"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Set

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

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

CS_TRIVIAL: float = 2.5   # ≤2.5: trivial, is_terminal=true
CS_MODERATE: float = 6.0  # ≤6.0: moderate; >6.0: complex

# Common protecting-group SMARTS for fg_density dimension.
# Kept local to avoid circular dependency on fg_detect (M2).
_PG_SMARTS: List[str] = [
    "[Si](C)(C)C",          # TMS / TBS / TIPS silyl ethers
    "C(=O)OC(C)(C)C",       # Boc
    "C(=O)OCc1ccccc1",      # Cbz
    "C(=O)C",               # Ac  (acetyl)
    "[CH2]c1ccccc1",        # Bn  (benzyl)
    "S(=O)(=O)c1ccc(C)cc1", # Ts  (tosyl)
    "C(=O)OC",              # methyl ester / Moc
    "C(OC)(OC)",            # acetal / ketal
]


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


def _dim_fg_density(mol: Chem.Mol) -> float:
    """FG density: count complexity-adding functional groups. Cap 2.5.

    Counts both protecting groups and synthetic-step-adding FGs
    (amides, sulfonamides, carbamates, ureas, etc.) that each typically
    require a dedicated coupling/formation step.
    """
    # Protecting groups (higher weight — they imply protect+deprotect steps)
    n_pg = 0
    for smarts_str in _PG_SMARTS:
        matches = smarts_match(mol, smarts_str)
        n_pg += len(matches)

    # Complexity-adding functional groups (each ≈ 1 synthetic step)
    _COMPLEX_FG = [
        "[NX3][CX3](=[OX1])[#6]",     # amide bond
        "[#6]S(=O)(=O)[NX3]",          # sulfonamide
        "[OX2][CX3](=[OX1])[NX3]",     # carbamate
        "[NX3][CX3](=[OX1])[NX3]",     # urea
        "[#6][CX3](=[OX1])[OX2][#6]",  # ester
        "[NX3][CX3](=[OX1])[OX2]",     # mixed: urethane / carbamate
        "[#6]C(=O)Cl",                  # acid chloride (reactive)
        "[#6]B([OH])[OH]",             # boronic acid
    ]
    n_fg = 0
    for smarts_str in _COMPLEX_FG:
        matches = smarts_match(mol, smarts_str)
        n_fg += len(matches)

    score = n_pg * 0.35 + n_fg * 0.25
    return min(score, 2.5)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_cs_score(smiles: str) -> Dict[str, Any]:
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

    # Compute raw dimension values
    raw_size = _dim_size(mol)
    raw_ring = _dim_ring(mol)
    raw_stereo = _dim_stereo(mol)
    raw_hetero = _dim_hetero(mol)
    raw_symmetry = _dim_symmetry(mol)   # positive magnitude
    raw_fg = _dim_fg_density(mol)

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
        "breakdown": {
            "size": round(c_size, 4),
            "ring": round(c_ring, 4),
            "stereo": round(c_stereo, 4),
            "hetero": round(c_hetero, 4),
            "symmetry": round(c_symmetry, 4),
            "fg_density": round(c_fg, 4),
        },
    }


def classify_complexity(smiles: str) -> Dict[str, Any]:
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
    result = compute_cs_score(smiles)

    # Propagate error
    if "ok" in result and result["ok"] is False:
        return result

    return {
        "classification": result["classification"],
        "is_terminal": result["is_terminal"],
        "is_heuristic": True,
    }


# ---------------------------------------------------------------------------
# Common FG SMARTS pairs for FG conversion bonus detection.
# Each tuple: (source_smarts, target_smarts, description)
# If intermediate has source FG and target has target FG, it suggests a
# plausible functional-group conversion step → small bonus.
# ---------------------------------------------------------------------------

_FG_CONVERSION_PAIRS: List[tuple] = [
    ("[OH]", "[OX2]C=O", "alcohol → ester"),
    ("[NH2]", "[NX3]C=O", "amine → amide"),
    ("C=C", "C-C", "alkene → alkane (reduction)"),
    ("[CH]=O", "C(=O)[OH]", "aldehyde → carboxylic acid"),
    ("[CH]=O", "C([OH])", "aldehyde → alcohol (reduction)"),
    ("C#N", "C(=O)[OH]", "nitrile → carboxylic acid"),
    ("C(=O)Cl", "C(=O)[OH]", "acyl chloride → carboxylic acid"),
    ("[Br]", "[OH]", "bromide → alcohol (substitution)"),
    ("[Cl]", "[OH]", "chloride → alcohol (substitution)"),
]


# ---------------------------------------------------------------------------
# score_progress helpers
# ---------------------------------------------------------------------------

def _count_pg_matches(mol: Chem.Mol) -> int:
    """Count total protecting-group matches on *mol* using _PG_SMARTS."""
    n = 0
    for smarts_str in _PG_SMARTS:
        n += len(smarts_match(mol, smarts_str))
    return n


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


def _fg_conversion_bonus(mol_inter: Chem.Mol, mol_target: Chem.Mol) -> float:
    """Return +0.02 if any FG conversion pair matches (intermediate has source,
    target has dest), else 0.0."""
    for src_smarts, dst_smarts, _desc in _FG_CONVERSION_PAIRS:
        if smarts_match(mol_inter, src_smarts) and smarts_match(mol_target, dst_smarts):
            return 0.02
    return 0.0


# ---------------------------------------------------------------------------
# Public API — progress scoring
# ---------------------------------------------------------------------------

def score_progress(intermediate: str, target: str) -> Dict[str, Any]:
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

    # 1. Tanimoto base
    tan_score = tanimoto(mol_inter, mol_target)

    # 2. Substructure bonus
    sub_bonus = 0.05 if _has_substructure_relation(mol_inter, mol_target) else 0.0

    # 3. Deprotection bonus: +0.03 per PG that intermediate has but target doesn't
    pg_inter = _count_pg_matches(mol_inter)
    pg_target = _count_pg_matches(mol_target)
    extra_pg = max(pg_inter - pg_target, 0)
    deprot_bonus = extra_pg * 0.03

    # 4. FG conversion bonus
    fg_bonus = _fg_conversion_bonus(mol_inter, mol_target)

    progress = min(1.0, tan_score + sub_bonus + deprot_bonus + fg_bonus)

    return {
        "progress_score": round(progress, 6),
        "tanimoto": round(tan_score, 6),
        "substructure_bonus": round(sub_bonus, 6),
        "deprotection_bonus": round(deprot_bonus, 6),
        "fg_conversion_bonus": round(fg_bonus, 6),
        "is_heuristic": True,
    }


def batch_score_progress(
    intermediates: List[str], target: str
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

    # Pre-compute target-side data
    fp_target = AllChem.GetMorganFingerprintAsBitVect(mol_target, radius=2, nBits=2048)
    pg_target = _count_pg_matches(mol_target)

    # Pre-compile FG conversion destination patterns for target
    _fg_dst_hits: List[bool] = []
    for _src, dst_smarts, _desc in _FG_CONVERSION_PAIRS:
        _fg_dst_hits.append(bool(smarts_match(mol_target, dst_smarts)))

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
        pg_inter = _count_pg_matches(mol_inter)
        extra_pg = max(pg_inter - pg_target, 0)
        deprot_bonus = extra_pg * 0.03

        # 4. FG conversion bonus (use pre-computed target dst hits)
        fg_bonus = 0.0
        for idx, (src_smarts, _dst, _desc) in enumerate(_FG_CONVERSION_PAIRS):
            if _fg_dst_hits[idx] and smarts_match(mol_inter, src_smarts):
                fg_bonus = 0.02
                break

        progress = min(1.0, tan_score + sub_bonus + deprot_bonus + fg_bonus)

        results.append({
            "progress_score": round(progress, 6),
            "tanimoto": round(tan_score, 6),
            "substructure_bonus": round(sub_bonus, 6),
            "deprotection_bonus": round(deprot_bonus, 6),
            "fg_conversion_bonus": round(fg_bonus, 6),
            "is_heuristic": True,
        })

    return results
