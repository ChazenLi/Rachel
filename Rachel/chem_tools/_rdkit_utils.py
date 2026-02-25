"""M0: RDKit 公共工具 — SMILES解析、SMARTS匹配、缓存、模板加载、常量定义。"""

from __future__ import annotations

import json
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# ---------------------------------------------------------------------------
# Constants for SMILES validation
# ---------------------------------------------------------------------------

_MAX_SMILES_LENGTH = 500
_MAX_NESTING_DEPTH = 20
_VALID_SMILES_CHARS = re.compile(r'^[A-Za-z0-9@+\-\[\]\(\)\.\#\=\:\%\/\\]+$')


# ---------------------------------------------------------------------------
# Core parsing & caching
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1024)
def parse_mol(smiles: str) -> Optional[Chem.Mol]:
    """SMILES → RDKit Mol with Sanitize. Returns None on failure.

    Results are cached via ``@lru_cache`` so repeated calls with the same
    SMILES string return the *same* object (identity-equal).
    """
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def canonical(smiles: str) -> Optional[str]:
    """Return the canonical SMILES string, or None if *smiles* is invalid."""
    mol = parse_mol(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def validate_smiles(smiles: str) -> Tuple[bool, str]:
    """Validate a SMILES string.

    Checks performed (in order):
    1. Non-empty after stripping whitespace
    2. Length ≤ 500
    3. Only legal SMILES characters
    4. Bracket / parenthesis nesting depth ≤ 20
    5. RDKit can parse it

    Returns
    -------
    (is_valid, reason) : tuple[bool, str]
        ``reason`` is ``""`` when valid, otherwise a human-readable message.
    """
    if not smiles or not smiles.strip():
        return False, "SMILES string is empty"

    if len(smiles) > _MAX_SMILES_LENGTH:
        return False, f"SMILES exceeds max length ({len(smiles)} > {_MAX_SMILES_LENGTH})"

    if not _VALID_SMILES_CHARS.match(smiles):
        return False, "SMILES contains invalid characters"

    # Check bracket / parenthesis nesting depth
    depth = 0
    max_depth = 0
    for ch in smiles:
        if ch in ("(", "["):
            depth += 1
            max_depth = max(max_depth, depth)
        elif ch in (")", "]"):
            depth -= 1

    if max_depth > _MAX_NESTING_DEPTH:
        return (
            False,
            f"SMILES nesting too deep ({max_depth} > {_MAX_NESTING_DEPTH})",
        )

    # Final check: RDKit parsing
    mol = parse_mol(smiles)
    if mol is None:
        return False, "RDKit cannot parse this SMILES"

    return True, ""


# ---------------------------------------------------------------------------
# SMARTS matching (with compiled-pattern cache)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=512)
def _compile_smarts(smarts: str) -> Optional[Chem.Mol]:
    """Compile a SMARTS string into an RDKit query Mol. Cached via LRU."""
    if not smarts:
        return None
    try:
        return Chem.MolFromSmarts(smarts)
    except Exception:
        return None


def smarts_match(mol: Chem.Mol, smarts: str) -> Tuple[Tuple[int, ...], ...]:
    """Return all substructure matches of *smarts* in *mol*.

    Uses an LRU-cached SMARTS compilation so repeated patterns are fast.
    Returns the same tuple-of-tuples as ``mol.GetSubstructMatches(pattern)``.
    """
    pattern = _compile_smarts(smarts)
    if pattern is None:
        return ()
    try:
        return mol.GetSubstructMatches(pattern)
    except Exception:
        return ()


# ---------------------------------------------------------------------------
# Tanimoto similarity
# ---------------------------------------------------------------------------

def tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    """Morgan fingerprint (radius=2, 2048-bit) Tanimoto similarity.

    Returns a float in [0.0, 1.0].
    """
    fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius=2, nBits=2048)
    fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius=2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


# ---------------------------------------------------------------------------
# Atom counting
# ---------------------------------------------------------------------------

def mol_formula_counter(mol: Chem.Mol) -> Counter:
    """Atom count including implicit hydrogens.

    Adds explicit Hs via ``Chem.AddHs``, then iterates every atom to build
    a ``Counter({'C': 6, 'H': 12, ...})``.
    """
    mol_h = Chem.AddHs(mol)
    cnt: Counter = Counter()
    for atom in mol_h.GetAtoms():
        cnt[atom.GetSymbol()] += 1
    return cnt


# ---------------------------------------------------------------------------
# Template loading
# ---------------------------------------------------------------------------

_TEMPLATES_DIR = Path(__file__).parent / "templates"


def _strip_comments(obj):
    """Recursively remove all ``__comment`` keys from a JSON-decoded object."""
    if isinstance(obj, dict):
        return {k: _strip_comments(v) for k, v in obj.items() if k != "__comment"}
    if isinstance(obj, list):
        return [_strip_comments(item) for item in obj]
    return obj


def load_template(filename: str) -> dict:
    """Load a JSON template from the ``templates/`` directory.

    Uses ``pathlib`` to locate the file relative to *this* module.
    All keys named ``"__comment"`` (at any nesting level) are stripped from
    the returned dictionary.
    """
    path = _TEMPLATES_DIR / filename
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return _strip_comments(data)


# ---------------------------------------------------------------------------
# Constants (migrated from old-Rachel/chem_analyze/graph_constants.py)
# ---------------------------------------------------------------------------

ELECTRONEGATIVITY: Dict[str, float] = {
    "H": 2.20, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
    "P": 2.19, "S": 2.58, "Cl": 3.16, "Se": 2.55, "Br": 2.96,
    "I": 2.66, "B": 2.04, "Si": 1.90, "Zn": 1.65, "Mg": 1.31,
    "Sn": 1.96, "Li": 0.98, "Na": 0.93, "K": 0.82,
    "Cu": 1.90, "Pd": 2.20, "Ni": 1.91, "Fe": 1.83,
}

VALENCE_ELECTRONS: Dict[str, int] = {
    "H": 1, "C": 4, "N": 5, "O": 6, "F": 7,
    "P": 5, "S": 6, "Cl": 7, "Se": 6, "Br": 7,
    "I": 7, "B": 3, "Si": 4, "Zn": 2, "Mg": 2,
    "Sn": 4, "Li": 1, "Na": 1, "K": 1,
    "Cu": 1, "Pd": 0, "Ni": 2, "Fe": 2,
}

FRAGMENT_TEMPLATES: Dict[str, str] = {
    "H": "[H]",
    "Br": "[Br]",
    "Cl": "[Cl]",
    "I": "[I]",
    "F": "[F]",
    "B(O)O": "B(O)O",
    "[Zn]Cl": "[Zn]Cl",
    "[Mg]Br": "[Mg]Br",
    "[Sn](C)(C)C": "[Sn](C)(C)C",
    "[Li]": "[Li]",
    "=O": "[O]",
    "OH": "[OH]",
    "CHO": "[CH]=O",
    "COOH": "C(=O)[OH]",
    "COCl": "C(=O)Cl",
    "NH2": "[NH2]",
    "NHMe": "[NH]C",
    "SO2Cl": "S(=O)(=O)Cl",
    "SiMe3": "[Si](C)(C)C",
    "CN": "C#N",
    "SH": "[SH]",
    "OAc": "OC(=O)C",
}
