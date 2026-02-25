"""M3: 反应模板扫描与匹配 — 扫描reactions.json中的模板，识别可断键位点。"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from ._rdkit_utils import canonical, load_template, parse_mol, smarts_match

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class TemplateMatch:
    """反应模板匹配结果。"""

    template_id: str                          # e.g. "C08_suzuki_retro"
    name: str                                 # e.g. "Suzuki Coupling (Retro)"
    category: str                             # e.g. "cc_formation"
    type: str                                 # "retro" | "forward"
    prereq_smarts: str                        # prerequisite SMARTS
    rxn_smarts: str                           # reaction SMARTS
    matched_atoms: List[Tuple[int, ...]]      # prereq match atom indices
    confidence: float                         # [0, 1], is_heuristic
    broken_bonds: List[Tuple[int, int]] = field(default_factory=list)
    incompatible_groups: List[str] = field(default_factory=list)
    conditions: Optional[str] = None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def scan_applicable_reactions(
    smiles: str,
    mode: str = "retro",
) -> Dict[str, Any]:
    """Scan reactions.json templates and return matches for *smiles*.

    Parameters
    ----------
    smiles : str
        Target molecule SMILES.
    mode : str
        ``"retro"``, ``"forward"``, or ``"both"`` — filters templates by
        their ``type`` field.

    Returns
    -------
    dict
        ``{"ok": True, "matches": [...], "by_category": {...}, "summary": {...}}``
        or ``{"ok": False, "error": "...", "input": smiles}`` on invalid input.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    can = canonical(smiles)
    n_heavy = mol.GetNumHeavyAtoms()

    # Load all reaction templates
    try:
        templates = load_template("reactions.json")
    except Exception as exc:
        return {"ok": False, "error": f"failed to load reactions.json: {exc}", "input": smiles}

    matches: List[TemplateMatch] = []

    for tpl_id, tpl in templates.items():
        if not isinstance(tpl, dict):
            continue

        tpl_type = tpl.get("type", "")

        # Filter by mode
        if mode == "retro" and tpl_type != "retro":
            continue
        if mode == "forward" and tpl_type != "forward":
            continue
        # mode == "both" → accept all types

        prereq = tpl.get("prereq", "")
        if not prereq:
            continue

        # Match prereq SMARTS against molecule
        matched = smarts_match(mol, prereq)
        if not matched:
            continue

        # Compute confidence = matched_atoms_count / heavy_atoms
        unique_atoms = set()
        for match_tuple in matched:
            unique_atoms.update(match_tuple)
        confidence = len(unique_atoms) / n_heavy if n_heavy > 0 else 0.0
        confidence = round(min(confidence, 1.0), 4)

        rxn = tpl.get("rxn", "")
        incompat = tpl.get("incompatible_groups", [])
        conditions = tpl.get("conditions") or tpl.get("note")

        tm = TemplateMatch(
            template_id=tpl_id,
            name=tpl.get("name", tpl_id),
            category=tpl.get("category", ""),
            type=tpl_type,
            prereq_smarts=prereq,
            rxn_smarts=rxn,
            matched_atoms=list(matched),
            confidence=confidence,
            broken_bonds=[],          # populated by extract_broken_bonds (task 7.2)
            incompatible_groups=list(incompat) if isinstance(incompat, list) else [],
            conditions=conditions,
        )
        matches.append(tm)

    # Sort by confidence descending
    matches.sort(key=lambda m: -m.confidence)

    # Group by category
    by_category: Dict[str, List[TemplateMatch]] = {}
    for m in matches:
        by_category.setdefault(m.category, []).append(m)

    categories = sorted(by_category.keys())

    return {
        "ok": True,
        "smiles": can,
        "matches": matches,
        "by_category": by_category,
        "summary": {
            "total": len(matches),
            "categories": categories,
        },
    }


# ---------------------------------------------------------------------------
# Broken-bond extraction from rxn SMARTS
# ---------------------------------------------------------------------------

# Module-level cache: rxn_smarts -> frozenset of (map_i, map_j)
_broken_pairs_cache: Dict[str, FrozenSet[Tuple[int, int]]] = {}


def _extract_broken_map_pairs(rxn_smarts: str) -> FrozenSet[Tuple[int, int]]:
    """Extract atom-map-number pairs whose bond is broken in *rxn_smarts*.

    Compares reactant-side and product-side bond topology:
      - Bond present in reactant but absent in product → broken bond.
      - Two mapped atoms in the same reactant template but in different
        product templates → bond broken (molecule split).

    FGI-type templates (bond order change only) yield an empty set.
    Results are cached per *rxn_smarts* string.
    """
    if rxn_smarts in _broken_pairs_cache:
        return _broken_pairs_cache[rxn_smarts]

    result: Set[Tuple[int, int]] = set()

    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception:
        _broken_pairs_cache[rxn_smarts] = frozenset()
        return frozenset()

    if rxn is None:
        _broken_pairs_cache[rxn_smarts] = frozenset()
        return frozenset()

    reactant_templates = [
        rxn.GetReactantTemplate(i)
        for i in range(rxn.GetNumReactantTemplates())
    ]
    product_templates = [
        rxn.GetProductTemplate(i)
        for i in range(rxn.GetNumProductTemplates())
    ]

    def _get_mapped_bonds(templates: list) -> Set[Tuple[int, int]]:
        bonds: Set[Tuple[int, int]] = set()
        for tmpl in templates:
            map_lookup: Dict[int, int] = {}
            for atom in tmpl.GetAtoms():
                mn = atom.GetAtomMapNum()
                if mn > 0:
                    map_lookup[atom.GetIdx()] = mn
            for bond in tmpl.GetBonds():
                m1 = map_lookup.get(bond.GetBeginAtomIdx())
                m2 = map_lookup.get(bond.GetEndAtomIdx())
                if m1 and m2:
                    bonds.add((min(m1, m2), max(m1, m2)))
        return bonds

    def _get_map_to_template_idx(templates: list) -> Dict[int, int]:
        m2t: Dict[int, int] = {}
        for tidx, tmpl in enumerate(templates):
            for atom in tmpl.GetAtoms():
                mn = atom.GetAtomMapNum()
                if mn > 0:
                    m2t[mn] = tidx
        return m2t

    reactant_bonds = _get_mapped_bonds(reactant_templates)
    product_bonds = _get_mapped_bonds(product_templates)

    # Strategy 1: bond in reactant but not in product
    for pair in reactant_bonds:
        if pair not in product_bonds:
            result.add(pair)

    # Strategy 2: same reactant template, different product templates
    reactant_m2t = _get_map_to_template_idx(reactant_templates)
    product_m2t = _get_map_to_template_idx(product_templates)
    shared_maps = set(reactant_m2t.keys()) & set(product_m2t.keys())

    for m1 in shared_maps:
        for m2 in shared_maps:
            if m1 >= m2:
                continue
            pair = (m1, m2)
            if pair in result:
                continue
            if (pair in reactant_bonds
                    and reactant_m2t.get(m1) == reactant_m2t.get(m2)
                    and product_m2t.get(m1) != product_m2t.get(m2)):
                result.add(pair)

    frozen = frozenset(result)
    _broken_pairs_cache[rxn_smarts] = frozen
    return frozen


def extract_broken_bonds(
    rxn_smarts: str,
    prereq_smarts: str,
    mol: Chem.Mol,
) -> Set[Tuple[int, int]]:
    """Extract broken bonds from *rxn_smarts* mapped to real atom indices in *mol*.

    Steps:
      1. Parse rxn SMARTS to find atom-map-number pairs whose bond is broken.
      2. Build a mapping from atom map numbers to prereq SMARTS atom indices.
      3. Use prereq substructure match to map prereq atom indices to real
         molecule atom indices.
      4. Verify each mapped pair actually has a bond in *mol*.

    Parameters
    ----------
    rxn_smarts : str
        Reaction SMARTS (``reactant>>product``).
    prereq_smarts : str
        Prerequisite SMARTS used for substructure matching.
    mol : Chem.Mol
        Target molecule.

    Returns
    -------
    Set[Tuple[int, int]]
        Set of ``(i, j)`` tuples (``i < j``) representing broken bonds in
        the target molecule.  Empty set on any parsing failure.
    """
    if not rxn_smarts or not prereq_smarts or mol is None:
        return set()

    try:
        broken_pairs = _extract_broken_map_pairs(rxn_smarts)
    except Exception:
        logger.warning("Failed to extract broken map pairs from rxn SMARTS: %s", rxn_smarts)
        return set()

    if not broken_pairs:
        return set()

    # Parse prereq SMARTS
    try:
        prereq_mol = Chem.MolFromSmarts(prereq_smarts)
    except Exception:
        logger.warning("Failed to parse prereq SMARTS: %s", prereq_smarts)
        return set()

    if prereq_mol is None:
        return set()

    # Build atom map number → prereq atom index
    map_to_prereq_idx: Dict[int, int] = {}
    for atom in prereq_mol.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn > 0:
            map_to_prereq_idx[mn] = atom.GetIdx()

    # Substructure match prereq against target molecule
    try:
        matches = mol.GetSubstructMatches(prereq_mol)
    except Exception:
        return set()

    if not matches:
        return set()

    real_bonds: Set[Tuple[int, int]] = set()
    for match in matches:
        for m1, m2 in broken_pairs:
            pi1 = map_to_prereq_idx.get(m1)
            pi2 = map_to_prereq_idx.get(m2)
            if pi1 is None or pi2 is None:
                continue
            if pi1 >= len(match) or pi2 >= len(match):
                continue
            ra1 = match[pi1]
            ra2 = match[pi2]
            # Verify bond exists in target molecule
            if mol.GetBondBetweenAtoms(ra1, ra2) is not None:
                real_bonds.add((min(ra1, ra2), max(ra1, ra2)))

    return real_bonds


# ---------------------------------------------------------------------------
# Bond priority helpers
# ---------------------------------------------------------------------------

_BOND_PRIORITY: Dict[str, float] = {
    "SINGLE": 1.0,
    "DOUBLE": 0.8,
    "TRIPLE": 0.6,
    "AROMATIC": 0.4,
}


def _get_bond_priority(mol: Chem.Mol, i: int, j: int) -> float:
    """Return heuristic bond priority based on bond type and ring membership."""
    bond = mol.GetBondBetweenAtoms(i, j)
    if bond is None:
        return 0.0
    bt = str(bond.GetBondType())
    # RDKit returns e.g. "Chem.rdchem.BondType.SINGLE" or "SINGLE"
    for key, val in _BOND_PRIORITY.items():
        if key in bt:
            priority = val
            break
    else:
        priority = 0.5
    # In-ring penalty
    if bond.IsInRing():
        priority *= 0.7
    return priority


# ---------------------------------------------------------------------------
# Strategic bond identification (no template match)
# ---------------------------------------------------------------------------

def _find_strategic_bonds_no_template(
    mol: Chem.Mol,
    matched_bond_set: Set[Tuple[int, int]],
) -> List[Dict[str, Any]]:
    """Identify strategic bonds that have no template match.

    Heuristic: single bonds between two ring atoms (inter-ring linkers),
    or single bonds connecting a ring atom to a non-ring heavy atom with
    degree >= 2.
    """
    n_heavy = mol.GetNumHeavyAtoms()
    strategic: List[Dict[str, Any]] = []

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        key = (min(a1, a2), max(a1, a2))

        if key in matched_bond_set:
            continue

        # Only single bonds
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        # Skip bonds inside rings
        if bond.IsInRing():
            continue

        atom1 = mol.GetAtomWithIdx(a1)
        atom2 = mol.GetAtomWithIdx(a2)

        # Heuristic: inter-ring linker or ring-to-substituent
        is_strategic = False
        if atom1.IsInRing() and atom2.IsInRing():
            is_strategic = True  # linker between two ring systems
        elif (atom1.IsInRing() or atom2.IsInRing()):
            # ring-to-substituent: only if substituent has degree >= 2
            non_ring = atom2 if atom1.IsInRing() else atom1
            if non_ring.GetDegree() >= 2:
                is_strategic = True

        if is_strategic:
            strategic.append({
                "atoms": list(key),
                "bond_type": "SINGLE",
                "in_ring": False,
                "note": "strategic bond without template match",
            })

    return strategic


# ---------------------------------------------------------------------------
# find_disconnectable_bonds — public API
# ---------------------------------------------------------------------------

def find_disconnectable_bonds(smiles: str) -> Dict[str, Any]:
    """Identify disconnectable bonds with associated retro-synthetic templates.

    Workflow:
      1. Call ``scan_applicable_reactions(smiles, mode="retro")``.
      2. For each match, call ``extract_broken_bonds`` to get broken bonds.
      3. Group templates by bond ``(i, j)``.
      4. For each bond compute ``heuristic_score = sqrt(bond_priority × best_confidence)``.
      5. Identify strategic bonds without template matches.

    Parameters
    ----------
    smiles : str
        Target molecule SMILES.

    Returns
    -------
    dict
        ``{"ok": True, "bonds": [...], "unmatched_bonds": [...]}``
        or ``{"ok": False, "error": "...", "input": smiles}`` on failure.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    can = canonical(smiles)
    n_atoms = mol.GetNumAtoms()

    # Step 1: scan retro templates
    scan_result = scan_applicable_reactions(smiles, mode="retro")
    if not scan_result.get("ok"):
        return scan_result

    matches: List[TemplateMatch] = scan_result.get("matches", [])

    # Step 2 & 3: extract broken bonds and group by bond
    bond_templates: Dict[Tuple[int, int], List[Dict[str, Any]]] = {}

    for tm in matches:
        broken = extract_broken_bonds(tm.rxn_smarts, tm.prereq_smarts, mol)
        # Also update the TemplateMatch's broken_bonds field
        tm.broken_bonds = sorted(broken)

        for bond_key in broken:
            entry = {
                "template_id": tm.template_id,
                "name": tm.name,
                "category": tm.category,
                "confidence": tm.confidence,
                "rxn_smarts": tm.rxn_smarts,
                "conditions": tm.conditions,
            }
            bond_templates.setdefault(bond_key, []).append(entry)

    # Step 4: compute heuristic_score per bond
    bonds_out: List[Dict[str, Any]] = []
    matched_bond_set: Set[Tuple[int, int]] = set()

    for bond_key, templates in bond_templates.items():
        i, j = bond_key
        matched_bond_set.add(bond_key)

        bp = _get_bond_priority(mol, i, j)
        best_conf = max(t["confidence"] for t in templates)
        score = math.sqrt(max(bp, 0.0) * max(best_conf, 0.0))
        score = round(score, 4)

        bond_obj = mol.GetBondBetweenAtoms(i, j)
        bt = "SINGLE"
        in_ring = False
        if bond_obj is not None:
            for key in _BOND_PRIORITY:
                if key in str(bond_obj.GetBondType()):
                    bt = key
                    break
            in_ring = bond_obj.IsInRing()

        bonds_out.append({
            "atoms": [i, j],
            "bond_type": bt,
            "in_ring": in_ring,
            "heuristic_score": score,
            "is_heuristic": True,
            "bond_priority": round(bp, 4),
            "best_template_confidence": round(best_conf, 4),
            "templates": templates,
        })

    # Sort by heuristic_score descending
    bonds_out.sort(key=lambda b: -b["heuristic_score"])

    # Step 5: unmatched strategic bonds
    unmatched = _find_strategic_bonds_no_template(mol, matched_bond_set)

    return {
        "ok": True,
        "smiles": can,
        "bonds": bonds_out,
        "unmatched_bonds": unmatched,
        "summary": {
            "total_bonds_with_templates": len(bonds_out),
            "total_unmatched_strategic": len(unmatched),
        },
    }
