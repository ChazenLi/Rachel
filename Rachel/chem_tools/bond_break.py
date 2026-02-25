"""M4: 断键执行与前体生成 — 接收断键决策并在RDKit Mol上执行操作，生成前体SMILES。"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, RWMol

from ._rdkit_utils import FRAGMENT_TEMPLATES, canonical, load_template, parse_mol
from .template_scan import scan_applicable_reactions

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class BreakResult:
    """断键执行结果。"""

    success: bool
    precursors: List[str] = field(default_factory=list)
    operation_log: List[str] = field(default_factory=list)
    error: Optional[str] = None
    is_fallback: bool = False


# ---------------------------------------------------------------------------
# Three-level SanitizeMol fallback
# ---------------------------------------------------------------------------

def _sanitize_mol_fallback(mol: Chem.Mol) -> Tuple[bool, Optional[str]]:
    """Apply three-level SanitizeMol fallback strategy.

    Returns (success, error_message).
    """
    # Level 1: Standard SanitizeMol
    try:
        Chem.SanitizeMol(mol)
        return True, None
    except Exception as e1:
        logger.debug("SanitizeMol Level 1 failed: %s", e1)

    # Level 2: SANITIZE_ALL ^ SANITIZE_PROPERTIES
    try:
        Chem.SanitizeMol(
            mol,
            Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
        )
        return True, None
    except Exception as e2:
        logger.debug("SanitizeMol Level 2 failed: %s", e2)

    # Level 3: Only Kekulize + SetAromaticity
    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
        Chem.SetAromaticity(mol)
        return True, None
    except Exception as e3:
        logger.debug("SanitizeMol Level 3 failed: %s", e3)
        return False, f"SanitizeMol failed at all 3 levels: {e3}"


# ---------------------------------------------------------------------------
# Helper: sanitize a product set from RunReactants
# ---------------------------------------------------------------------------

def _sanitize_product_set(products) -> Optional[List[str]]:
    """Convert RunReactants product Mols to canonical SMILES list.

    Applies the three-level sanitize fallback per product molecule.
    Returns None if any product cannot be sanitized.
    """
    smiles_list: List[str] = []
    for prod_mol in products:
        ok, _err = _sanitize_mol_fallback(prod_mol)
        if not ok:
            return None
        try:
            smi = Chem.MolToSmiles(prod_mol, canonical=True)
            if smi:
                smiles_list.append(smi)
            else:
                return None
        except Exception:
            return None
    return smiles_list if smiles_list else None


# ---------------------------------------------------------------------------
# Helper: find retro templates for a reaction_type via scan
# ---------------------------------------------------------------------------


def _find_templates_for_reaction(
    smiles: str,
    reaction_type: str,
    bond: Tuple[int, int],
) -> List[Dict[str, Any]]:
    """Find retro templates matching *reaction_type* that break *bond*.

    Uses scan_applicable_reactions to get all retro matches, then filters
    by reaction_type.

    Matching strategy (in order):
      1. Exact template_id match (case-insensitive).
      2. Normalized substring match — strip parenthetical suffixes like
         ``"(Retro, ester amide)"`` and common prefixes like ``"C08_"``,
         then check substring containment against template_id, name,
         and category.
      3. Token overlap — split both strings into word tokens and check
         if ≥50% of the shorter token set appears in the longer one.
    """
    scan = scan_applicable_reactions(smiles, mode="retro")
    if not scan.get("ok"):
        return []

    # --- Normalize the input reaction_type ---
    raw = reaction_type.strip()
    raw_lower = raw.lower()

    # Strip parenthetical suffix: "Amide Coupling EDC HATU (Retro, ester amide)" → "Amide Coupling EDC HATU"
    paren_idx = raw.find("(")
    if paren_idx > 0:
        raw = raw[:paren_idx].strip()

    key = raw.lower().replace("-", "_").replace(" ", "_")
    # Also strip trailing _retro for compat
    key_no_retro = re.sub(r"_retro$", "", key)

    # Token set for fallback matching
    key_tokens = set(re.split(r"[\s_\-,()]+", raw_lower))
    key_tokens.discard("")
    # Remove very common noise tokens
    key_tokens -= {"retro", "forward", "reaction", "the", "a", "an", "of"}

    candidates: List[Dict[str, Any]] = []

    for tm in scan.get("matches", []):
        tid_lower = tm.template_id.lower()
        name_lower = tm.name.lower()
        cat_lower = tm.category.lower()

        matched = False

        # Strategy 1: exact template_id match
        if raw_lower == tid_lower:
            matched = True

        # Strategy 2: normalized substring match
        if not matched:
            # Strip prefix like "G14_" from template_id for comparison
            tid_stem = re.sub(r"^[A-Z]\d+_", "", tm.template_id, flags=re.IGNORECASE).lower()
            tid_stem_no_retro = re.sub(r"_retro$", "", tid_stem)

            if (
                key in tid_lower
                or key_no_retro in tid_lower
                or tid_stem_no_retro in key
                or key_no_retro in tid_stem_no_retro
                or tid_stem_no_retro in key_no_retro
                or key in name_lower
                or key_no_retro in name_lower
                or key in cat_lower
            ):
                matched = True

        # Strategy 3: token overlap (≥50% of shorter set)
        if not matched and key_tokens:
            tm_tokens = set(re.split(r"[\s_\-,()]+", f"{tid_lower} {name_lower} {cat_lower}"))
            tm_tokens.discard("")
            tm_tokens -= {"retro", "forward", "reaction", "the", "a", "an", "of"}

            overlap = key_tokens & tm_tokens
            shorter_len = min(len(key_tokens), len(tm_tokens))
            if shorter_len > 0 and len(overlap) / shorter_len >= 0.5:
                matched = True

        if matched:
            candidates.append({
                "template_id": tm.template_id,
                "rxn_smarts": tm.rxn_smarts,
                "prereq_smarts": tm.prereq_smarts,
                "name": tm.name,
                "category": tm.category,
            })

    return candidates




# ---------------------------------------------------------------------------
# Priority 2: rxn SMARTS template execution
# ---------------------------------------------------------------------------

def _try_template_execution(
    mol: Chem.Mol,
    bond: Tuple[int, int],
    templates: List[Dict[str, Any]],
    op_log: List[str],
) -> Optional[BreakResult]:
    """Try executing rxn SMARTS templates on the specified bond.

    For each template:
      1. Parse rxn SMARTS via AllChem.ReactionFromSmarts
      2. RunReactants on the molecule
      3. Filter product sets to find one that breaks the target bond
      4. Sanitize and return precursors

    Returns BreakResult on success, None if no template works.
    """
    target = (min(bond[0], bond[1]), max(bond[0], bond[1]))

    for tpl in templates:
        rxn_smarts = tpl.get("rxn_smarts", "")
        prereq_smarts = tpl.get("prereq_smarts", "")
        tpl_name = tpl.get("name", tpl.get("template_id", "unknown"))

        if not rxn_smarts:
            continue

        try:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        except Exception:
            op_log.append(f"Template {tpl_name}: rxn SMARTS parse failed")
            continue

        if rxn is None:
            op_log.append(f"Template {tpl_name}: rxn SMARTS returned None")
            continue

        try:
            product_sets = rxn.RunReactants((mol,))
        except Exception as exc:
            op_log.append(f"Template {tpl_name}: RunReactants failed: {exc}")
            continue

        if not product_sets:
            op_log.append(f"Template {tpl_name}: no products generated")
            continue

        # Try each product set — check if it breaks the target bond
        for prods in product_sets:
            # Use react_atom_idx to verify the correct bond was broken
            mapno_to_react_idx: Dict[int, int] = {}
            for prod_mol in prods:
                for atom in prod_mol.GetAtoms():
                    if atom.HasProp("old_mapno") and atom.HasProp("react_atom_idx"):
                        mapno_to_react_idx[atom.GetIntProp("old_mapno")] = (
                            atom.GetIntProp("react_atom_idx")
                        )

            # If we have mapping info, verify the bond
            if mapno_to_react_idx:
                # Extract broken map pairs from prereq
                from .template_scan import _extract_broken_map_pairs

                broken_pairs = _extract_broken_map_pairs(rxn_smarts)
                is_target = False
                for m1, m2 in broken_pairs:
                    ri1 = mapno_to_react_idx.get(m1)
                    ri2 = mapno_to_react_idx.get(m2)
                    if ri1 is not None and ri2 is not None:
                        if (min(ri1, ri2), max(ri1, ri2)) == target:
                            is_target = True
                            break
                if not is_target:
                    continue

            precursors = _sanitize_product_set(prods)
            if precursors and len(precursors) >= 2:
                op_log.append(
                    f"Template {tpl_name} succeeded: {precursors}"
                )
                return BreakResult(
                    success=True,
                    precursors=precursors,
                    operation_log=op_log,
                )

        op_log.append(f"Template {tpl_name}: no valid product set for bond {target}")

    return None


# ---------------------------------------------------------------------------
# Priority 3: Simple break + add H fallback
# ---------------------------------------------------------------------------

def _simple_break(
    mol: Chem.Mol,
    bond: Tuple[int, int],
    op_log: List[str],
) -> BreakResult:
    """Simple bond break: RemoveBond → GetMolFrags → add H to broken atoms → SanitizeMol.

    This is the fallback when no template matches.
    """
    i, j = bond
    rw = RWMol(mol)

    # Remove the bond
    rw.RemoveBond(i, j)
    op_log.append(f"Simple break: removed bond {i}-{j}")

    # Add explicit H to the broken atoms to satisfy valence
    for atom_idx in (i, j):
        atom = rw.GetAtomWithIdx(atom_idx)
        current_hs = atom.GetNumExplicitHs()
        atom.SetNumExplicitHs(current_hs + 1)
        atom.SetNoImplicit(False)

    op_log.append(f"Added H to atoms {i} and {j}")

    # Get the modified mol
    try:
        edited_mol = rw.GetMol()
    except Exception as exc:
        return BreakResult(
            success=False,
            operation_log=op_log,
            error=f"Failed to get mol after bond removal: {exc}",
        )

    # Three-level sanitize fallback
    ok, err = _sanitize_mol_fallback(edited_mol)
    if not ok:
        return BreakResult(
            success=False,
            operation_log=op_log,
            error=err,
        )

    # Split into fragments
    try:
        frag_mols = Chem.GetMolFrags(edited_mol, asMols=True, sanitizeFrags=True)
    except Exception:
        try:
            frag_mols = Chem.GetMolFrags(edited_mol, asMols=True, sanitizeFrags=False)
        except Exception as exc:
            return BreakResult(
                success=False,
                operation_log=op_log,
                error=f"GetMolFrags failed: {exc}",
            )

    if not frag_mols:
        return BreakResult(
            success=False,
            operation_log=op_log,
            error="No fragments generated after bond removal",
        )

    # Convert fragments to SMILES
    precursors: List[str] = []
    for frag in frag_mols:
        # Sanitize each fragment individually
        ok_frag, err_frag = _sanitize_mol_fallback(frag)
        if not ok_frag:
            return BreakResult(
                success=False,
                operation_log=op_log,
                error=f"Fragment sanitize failed: {err_frag}",
            )
        try:
            smi = Chem.MolToSmiles(frag, canonical=True)
        except Exception as exc:
            return BreakResult(
                success=False,
                operation_log=op_log,
                error=f"Fragment to SMILES failed: {exc}",
            )
        if not smi:
            return BreakResult(
                success=False,
                operation_log=op_log,
                error="Fragment produced empty SMILES",
            )
        # Verify round-trip
        check = Chem.MolFromSmiles(smi)
        if check is None:
            return BreakResult(
                success=False,
                operation_log=op_log,
                error=f"Generated SMILES invalid: {smi}",
            )
        precursors.append(smi)

    op_log.append(f"Simple break produced: {precursors}")
    return BreakResult(
        success=True,
        precursors=precursors,
        operation_log=op_log,
        is_fallback=True,
    )


# ---------------------------------------------------------------------------
# Public API: execute_disconnection
# ---------------------------------------------------------------------------

def execute_disconnection(
    smiles: str,
    bond: Tuple[int, int],
    reaction_type: str,
    custom_fragments: Optional[List[Tuple[str, str]]] = None,
) -> BreakResult:
    """Execute a single bond disconnection on *smiles* at *bond*.

    Three-level priority:
      1. **custom_fragments** — If provided, use LLM-specified fragment SMILES
         directly (validate each SMILES).
      2. **rxn SMARTS template** — Find matching retro templates for
         *reaction_type* and execute via ``AllChem.ReactionFromSmarts``.
      3. **Simple break + add H** — Fallback: ``RWMol.RemoveBond`` →
         ``GetMolFrags`` → add H to broken atoms → ``SanitizeMol``.

    Parameters
    ----------
    smiles : str
        Target molecule SMILES.
    bond : Tuple[int, int]
        Atom index pair ``(i, j)`` of the bond to break.
    reaction_type : str
        Reaction type name (e.g. ``"suzuki"``, ``"amide_formation"``).
    custom_fragments : Optional[List[Tuple[str, str]]]
        LLM-specified fragment SMILES pairs. Each tuple is
        ``(fragment_smiles_1, fragment_smiles_2)``.

    Returns
    -------
    BreakResult
        Contains ``success``, ``precursors``, ``operation_log``, ``error``.
    """
    op_log: List[str] = []

    # --- Validate input SMILES ---
    mol = parse_mol(smiles)
    if mol is None:
        return BreakResult(
            success=False,
            operation_log=op_log,
            error=f"invalid SMILES: {smiles}",
        )

    can = canonical(smiles)
    op_log.append(f"Target: {can}")

    # --- Validate bond ---
    i, j = int(bond[0]), int(bond[1])
    n_atoms = mol.GetNumAtoms()
    if i < 0 or j < 0 or i >= n_atoms or j >= n_atoms:
        return BreakResult(
            success=False,
            operation_log=op_log,
            error=f"atom index out of range: ({i}, {j}), molecule has {n_atoms} atoms",
        )

    if mol.GetBondBetweenAtoms(i, j) is None:
        return BreakResult(
            success=False,
            operation_log=op_log,
            error=f"no bond between atoms {i} and {j}",
        )

    op_log.append(f"Breaking bond: {i}-{j}, reaction_type={reaction_type}")

    # ===================================================================
    # Priority 1: custom_fragments
    # ===================================================================
    if custom_fragments is not None:
        op_log.append("Priority 1: using custom_fragments")
        all_precursors: List[str] = []
        for frag_pair in custom_fragments:
            for frag_smi in frag_pair:
                frag_mol = parse_mol(frag_smi)
                if frag_mol is None:
                    return BreakResult(
                        success=False,
                        operation_log=op_log,
                        error=f"invalid custom fragment SMILES: {frag_smi}",
                    )
                csmi = canonical(frag_smi)
                if csmi:
                    all_precursors.append(csmi)
                else:
                    all_precursors.append(frag_smi)

        if all_precursors:
            op_log.append(f"Custom fragments validated: {all_precursors}")
            return BreakResult(
                success=True,
                precursors=all_precursors,
                operation_log=op_log,
            )
        else:
            return BreakResult(
                success=False,
                operation_log=op_log,
                error="custom_fragments provided but no valid fragments found",
            )

    # ===================================================================
    # Priority 2: rxn SMARTS template execution
    # ===================================================================
    op_log.append("Priority 2: trying rxn SMARTS templates")
    templates = _find_templates_for_reaction(smiles, reaction_type, (i, j))

    if templates:
        op_log.append(f"Found {len(templates)} candidate templates for '{reaction_type}'")
        result = _try_template_execution(mol, (i, j), templates, op_log)
        if result is not None:
            return result
        op_log.append("All templates failed, falling back to simple break")
    else:
        op_log.append(f"No templates found for '{reaction_type}', falling back to simple break")

    # ===================================================================
    # Priority 3: Simple break + add H
    # ===================================================================
    op_log.append("Priority 3: simple break + add H fallback")
    return _simple_break(mol, (i, j), op_log)

def preview_disconnections(
    smiles: str,
    bond: Tuple[int, int],
) -> Dict[str, Any]:
    """Preview all template-based disconnection alternatives for a single bond.

    For each retro template that matches the target bond, execute the
    disconnection and collect the resulting precursors.  Returns a
    compact summary so the LLM can compare Suzuki vs Negishi vs Kumada
    (etc.) in one call without multiple round-trips.

    Parameters
    ----------
    smiles : str
        Target molecule SMILES.
    bond : Tuple[int, int]
        Atom index pair ``(i, j)`` of the bond to preview.

    Returns
    -------
    dict
        ``{"ok": True, "bond": [i, j], "alternatives": [...]}``

        Each alternative::

            {
                "template_id": str,
                "name": str,
                "category": str,
                "precursors": [str, ...],
                "leaving_groups": [str, ...],   # non-scaffold fragments
                "note": str,                    # template note for LLM context
                "incompatible_groups": [str],   # FG incompatibilities
            }

        Alternatives are sorted by template confidence (descending).
        On failure: ``{"ok": False, "error": "..."}``.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": f"invalid SMILES: {smiles}"}

    i, j = int(bond[0]), int(bond[1])
    n_atoms = mol.GetNumAtoms()
    if i < 0 or j < 0 or i >= n_atoms or j >= n_atoms:
        return {"ok": False, "error": f"atom index out of range: ({i}, {j})"}
    if mol.GetBondBetweenAtoms(i, j) is None:
        return {"ok": False, "error": f"no bond between atoms {i} and {j}"}

    # Scan all retro templates matching this molecule
    scan = scan_applicable_reactions(smiles, mode="retro")
    if not scan.get("ok"):
        return {"ok": False, "error": "template scan failed"}

    from .template_scan import extract_broken_bonds

    target = (min(i, j), max(i, j))
    alternatives: List[Dict[str, Any]] = []
    seen_precursor_sets: set = set()

    templates_data = _get_reactions()

    for tm in scan.get("matches", []):
        # Check if this template breaks our target bond
        broken = extract_broken_bonds(tm.rxn_smarts, tm.prereq_smarts, mol)
        if target not in broken:
            continue

        # Try executing
        result = try_retro_template(mol, (i, j), tm.template_id)
        if result is None or not result.success or not result.precursors:
            continue

        # Deduplicate by precursor set
        prec_key = frozenset(result.precursors)
        if prec_key in seen_precursor_sets:
            continue
        seen_precursor_sets.add(prec_key)

        # Extract metadata from reactions.json
        tpl_data = templates_data.get(tm.template_id, {})

        alternatives.append({
            "template_id": tm.template_id,
            "name": tm.name,
            "category": tm.category,
            "confidence": tm.confidence,
            "precursors": result.precursors,
            "note": tpl_data.get("note", ""),
            "incompatible_groups": tpl_data.get("incompatible_groups", []),
        })

    # Sort by confidence descending
    alternatives.sort(key=lambda a: -a.get("confidence", 0))

    return {
        "ok": True,
        "smiles": canonical(smiles),
        "bond": [i, j],
        "n_alternatives": len(alternatives),
        "alternatives": alternatives,
    }



# ---------------------------------------------------------------------------
# Module-level caches for resolve_template_ids
# ---------------------------------------------------------------------------

_REACTIONS_CACHE: Optional[Dict[str, Any]] = None
_NAME_INDEX_CACHE: Optional[Dict[str, List[str]]] = None


def _get_reactions() -> Dict[str, Any]:
    """Lazy-load reactions.json."""
    global _REACTIONS_CACHE
    if _REACTIONS_CACHE is not None:
        return _REACTIONS_CACHE
    try:
        _REACTIONS_CACHE = load_template("reactions.json")
    except Exception:
        _REACTIONS_CACHE = {}
    return _REACTIONS_CACHE


def _build_name_index() -> Dict[str, List[str]]:
    """Build a keyword → template_id list index from reactions.json."""
    global _NAME_INDEX_CACHE
    if _NAME_INDEX_CACHE is not None:
        return _NAME_INDEX_CACHE

    templates = _get_reactions()
    index: Dict[str, List[str]] = {}

    for tid, info in templates.items():
        if not isinstance(info, dict):
            continue
        if info.get("type") != "retro":
            continue

        # Extract stem from template_id: "C08_suzuki_retro" → "suzuki"
        stem = re.sub(r"^[A-Z]\d+_", "", tid)
        stem = re.sub(r"_retro$", "", stem)
        key = stem.lower()
        if key:
            index.setdefault(key, []).append(tid)

        # Use first word of name as alias
        name = info.get("name", "")
        if name:
            first_word = name.split()[0].lower().rstrip("(,")
            if first_word and first_word != key:
                index.setdefault(first_word, []).append(tid)

    _NAME_INDEX_CACHE = index
    return _NAME_INDEX_CACHE


# ---------------------------------------------------------------------------
# Public API: resolve_template_ids
# ---------------------------------------------------------------------------

def resolve_template_ids(reaction_name: str) -> List[str]:
    """Resolve a reaction name to a list of matching retro template IDs.

    Normalizes *reaction_name* (lowercase, strip spaces, replace hyphens
    with underscores) then searches reactions.json:

      1. Exact template_id match
      2. Keyword index lookup (stem extracted from template_id + first word
         of name field)
      3. Substring fuzzy match across all index keys

    Parameters
    ----------
    reaction_name : str
        Reaction name string (e.g. ``"suzuki"``, ``"Wittig"``,
        ``"C08_suzuki_retro"``).

    Returns
    -------
    List[str]
        Matching template IDs. Empty list if nothing matches.
    """
    templates = _get_reactions()

    # 1. Exact template_id
    if reaction_name in templates:
        info = templates[reaction_name]
        if isinstance(info, dict) and info.get("type") == "retro":
            return [reaction_name]

    # 2. Keyword index lookup
    index = _build_name_index()
    key = reaction_name.lower().strip().replace("-", "_").replace(" ", "_")
    # Strip trailing _retro for compat
    key = re.sub(r"_retro$", "", key)

    if key in index:
        return list(index[key])

    # 3. Substring fuzzy match
    matches: List[str] = []
    seen: Set[str] = set()
    for idx_key, tids in index.items():
        if key in idx_key or idx_key in key:
            for tid in tids:
                if tid not in seen:
                    seen.add(tid)
                    matches.append(tid)
    return matches


# ---------------------------------------------------------------------------
# Public API: try_retro_template
# ---------------------------------------------------------------------------

def try_retro_template(
    mol: Chem.Mol,
    bond: Tuple[int, int],
    template_id: str,
) -> Optional[BreakResult]:
    """Execute a specific retro template on *mol* at *bond*.

    Loads the template from reactions.json by *template_id*, verifies the
    prereq SMARTS matches and that the reaction centre covers the target
    bond, then runs ``AllChem.ReactionFromSmarts`` + ``RunReactants``.
    Only product sets that actually break the specified bond are kept.

    Parameters
    ----------
    mol : Chem.Mol
        Target molecule.
    bond : Tuple[int, int]
        Atom index pair ``(i, j)`` of the bond to break.
    template_id : str
        Template ID in reactions.json (e.g. ``"C08_suzuki_retro"``).

    Returns
    -------
    Optional[BreakResult]
        ``BreakResult`` on success, ``None`` if the template does not
        match or fails to produce valid products for the specified bond.
    """
    templates = _get_reactions()
    tpl = templates.get(template_id)
    if tpl is None or not isinstance(tpl, dict) or tpl.get("type") != "retro":
        return None

    prereq_smarts = tpl.get("prereq", "")
    rxn_smarts = tpl.get("rxn", "")
    if not prereq_smarts or not rxn_smarts:
        return None

    op_log: List[str] = [f"try_retro_template: {template_id}"]

    # --- prereq match ---
    prereq_pat = Chem.MolFromSmarts(prereq_smarts)
    if prereq_pat is None:
        op_log.append("prereq SMARTS parse failed")
        return None

    matches = mol.GetSubstructMatches(prereq_pat)
    if not matches:
        op_log.append("prereq SMARTS did not match")
        return None

    # --- extract broken map pairs ---
    from .template_scan import _extract_broken_map_pairs

    broken_map_pairs = _extract_broken_map_pairs(rxn_smarts)
    if not broken_map_pairs:
        op_log.append("no broken map pairs (FGI template)")
        return None

    # Build map_number → prereq internal atom idx
    map_to_pidx: Dict[int, int] = {}
    for atom in prereq_pat.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn > 0:
            map_to_pidx[mn] = atom.GetIdx()

    # --- verify target bond is a reaction centre ---
    target = (min(bond[0], bond[1]), max(bond[0], bond[1]))
    has_match = False
    for match in matches:
        for m1, m2 in broken_map_pairs:
            pi1 = map_to_pidx.get(m1)
            pi2 = map_to_pidx.get(m2)
            if pi1 is None or pi2 is None:
                continue
            if pi1 >= len(match) or pi2 >= len(match):
                continue
            ra1, ra2 = match[pi1], match[pi2]
            if (min(ra1, ra2), max(ra1, ra2)) == target:
                has_match = True
                break
        if has_match:
            break

    if not has_match:
        op_log.append(f"bond {target} is not a reaction centre for this template")
        return None

    # --- RunReactants + filter ---
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception:
        op_log.append("rxn SMARTS parse failed")
        return None

    if rxn is None:
        return None

    try:
        product_sets = rxn.RunReactants((mol,))
    except Exception as exc:
        op_log.append(f"RunReactants failed: {exc}")
        return None

    if not product_sets:
        op_log.append("no products generated")
        return None

    for prods in product_sets:
        # Use react_atom_idx to verify the correct bond was broken
        mapno_to_react_idx: Dict[int, int] = {}
        for prod_mol in prods:
            for atom in prod_mol.GetAtoms():
                if atom.HasProp("old_mapno") and atom.HasProp("react_atom_idx"):
                    mapno_to_react_idx[atom.GetIntProp("old_mapno")] = (
                        atom.GetIntProp("react_atom_idx")
                    )

        if mapno_to_react_idx:
            is_target = False
            for m1, m2 in broken_map_pairs:
                ri1 = mapno_to_react_idx.get(m1)
                ri2 = mapno_to_react_idx.get(m2)
                if ri1 is not None and ri2 is not None:
                    if (min(ri1, ri2), max(ri1, ri2)) == target:
                        is_target = True
                        break
            if not is_target:
                continue

        precursors = _sanitize_product_set(prods)
        if precursors and len(precursors) >= 2:
            op_log.append(f"template succeeded: {precursors}")
            return BreakResult(
                success=True,
                precursors=precursors,
                operation_log=op_log,
            )

    op_log.append("no valid product set for target bond")
    return None


# ---------------------------------------------------------------------------
# Public API: execute_fgi
# ---------------------------------------------------------------------------

def execute_fgi(
    smiles: str,
    template_id: str,
) -> BreakResult:
    """Execute a Functional Group Interconversion (FGI) retro template.

    FGI templates transform one functional group into another without
    breaking any bond (e.g. ketone ← alcohol, ArH ← ArBr).  They
    produce exactly **one** product (the precursor) rather than two
    fragments.

    Workflow:
      1. Load template from reactions.json by *template_id*.
      2. Verify prereq SMARTS matches the molecule.
      3. Run ``AllChem.ReactionFromSmarts`` + ``RunReactants``.
      4. Sanitize and return the single precursor.

    Parameters
    ----------
    smiles : str
        Target molecule SMILES.
    template_id : str
        Template ID in reactions.json (e.g. ``"A01_ketone_reduction_NaBH4_retro"``).
        Also accepts fuzzy names — resolved via ``resolve_template_ids``.

    Returns
    -------
    BreakResult
        ``success=True`` with a single-element ``precursors`` list on
        success.  ``is_fallback`` is always ``False`` for FGI results.
    """
    op_log: List[str] = [f"execute_fgi: {template_id}"]

    # --- Validate input ---
    mol = parse_mol(smiles)
    if mol is None:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"invalid SMILES: {smiles}",
        )

    can = canonical(smiles)
    op_log.append(f"Target: {can}")

    # --- Resolve template_id (support fuzzy names) ---
    templates = _get_reactions()
    tpl = templates.get(template_id)

    if tpl is None or not isinstance(tpl, dict) or tpl.get("type") != "retro":
        # Try fuzzy resolution
        resolved = resolve_template_ids(template_id)
        if not resolved:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"template not found: {template_id}",
            )
        template_id = resolved[0]
        tpl = templates.get(template_id)
        if tpl is None or not isinstance(tpl, dict):
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"resolved template invalid: {template_id}",
            )

    op_log.append(f"Using template: {template_id}")

    prereq_smarts = tpl.get("prereq", "")
    rxn_smarts = tpl.get("rxn", "")
    if not prereq_smarts or not rxn_smarts:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"template {template_id} missing prereq or rxn",
        )

    # --- Prereq match ---
    prereq_pat = Chem.MolFromSmarts(prereq_smarts)
    if prereq_pat is None:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"prereq SMARTS parse failed: {prereq_smarts}",
        )

    matches = mol.GetSubstructMatches(prereq_pat)
    if not matches:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"prereq SMARTS did not match molecule",
        )

    op_log.append(f"Prereq matched {len(matches)} site(s)")

    # --- Run reaction ---
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception as exc:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"rxn SMARTS parse failed: {exc}",
        )

    if rxn is None:
        return BreakResult(
            success=False, operation_log=op_log,
            error="rxn SMARTS returned None",
        )

    try:
        product_sets = rxn.RunReactants((mol,))
    except Exception as exc:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"RunReactants failed: {exc}",
        )

    if not product_sets:
        return BreakResult(
            success=False, operation_log=op_log,
            error="no products generated",
        )

    # --- Sanitize first valid product set ---
    for prods in product_sets:
        precursors = _sanitize_product_set(prods)
        if precursors:
            op_log.append(f"FGI succeeded: {precursors}")
            return BreakResult(
                success=True,
                precursors=precursors,
                operation_log=op_log,
            )

    return BreakResult(
        success=False, operation_log=op_log,
        error="all product sets failed sanitization",
    )


# ---------------------------------------------------------------------------
# Internal operation handlers for execute_operations
# ---------------------------------------------------------------------------

def _op_disconnect(
    rw: RWMol,
    params: Dict[str, Any],
    op_log: List[str],
) -> bool:
    """Handle 'disconnect' operation: remove a bond between two atoms."""
    bond_pair = params.get("bond") or params.get("atoms")
    if not bond_pair or len(bond_pair) != 2:
        op_log.append("disconnect: missing or invalid 'bond' param")
        return False

    i, j = int(bond_pair[0]), int(bond_pair[1])
    if rw.GetBondBetweenAtoms(i, j) is None:
        op_log.append(f"disconnect: no bond between {i} and {j}")
        return False

    rw.RemoveBond(i, j)

    # Add H to satisfy valence on broken atoms
    for idx in (i, j):
        atom = rw.GetAtomWithIdx(idx)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        atom.SetNoImplicit(False)

    op_log.append(f"disconnect: removed bond {i}-{j}")
    return True


def _op_add_fragment(
    rw: RWMol,
    params: Dict[str, Any],
    op_log: List[str],
) -> bool:
    """Handle 'add_fragment' operation: attach a fragment to an atom."""
    atom_idx = params.get("atom_idx")
    fragment = params.get("fragment") or params.get("group")
    if atom_idx is None or fragment is None:
        op_log.append("add_fragment: missing 'atom_idx' or 'fragment'")
        return False

    atom_idx = int(atom_idx)
    if atom_idx >= rw.GetNumAtoms():
        op_log.append(f"add_fragment: atom_idx {atom_idx} out of range")
        return False

    # Resolve fragment SMILES
    frag_smi = FRAGMENT_TEMPLATES.get(fragment, fragment)
    frag_mol = Chem.MolFromSmiles(frag_smi)
    if frag_mol is None:
        op_log.append(f"add_fragment: cannot parse fragment '{fragment}'")
        return False

    attach_idx = params.get("attach_idx", 0)

    # Combine molecules
    offset = rw.GetNumAtoms()
    combined = Chem.CombineMols(rw.GetMol(), frag_mol)
    rw_new = RWMol(combined)
    frag_attach = offset + int(attach_idx)

    bond_type_str = params.get("bond_type", "SINGLE")
    bond_type_map = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
    }
    bt = bond_type_map.get(bond_type_str.upper(), Chem.BondType.SINGLE)

    rw_new.AddBond(atom_idx, frag_attach, bt)

    # Reduce explicit H on the target atom if possible
    target_atom = rw_new.GetAtomWithIdx(atom_idx)
    hs = target_atom.GetNumExplicitHs()
    if hs > 0:
        target_atom.SetNumExplicitHs(hs - 1)

    # Copy back into rw (replace contents)
    # We need to return the new RWMol, so we modify in place via a trick:
    # Actually, we can't easily replace rw in-place. We'll use a different
    # approach — return the new mol and let the caller handle it.
    # For now, store on rw._replacement
    rw._replacement = rw_new.GetMol()  # type: ignore[attr-defined]
    op_log.append(f"add_fragment: attached '{fragment}' to atom {atom_idx}")
    return True


def _op_change_bond_order(
    rw: RWMol,
    params: Dict[str, Any],
    op_log: List[str],
) -> bool:
    """Handle 'change_bond_order' operation."""
    bond_pair = params.get("bond") or params.get("atoms")
    new_order = params.get("new_order")
    if not bond_pair or len(bond_pair) != 2 or new_order is None:
        op_log.append("change_bond_order: missing params")
        return False

    i, j = int(bond_pair[0]), int(bond_pair[1])
    bond = rw.GetBondBetweenAtoms(i, j)
    if bond is None:
        op_log.append(f"change_bond_order: no bond between {i} and {j}")
        return False

    order_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    new_bt = order_map.get(int(new_order))
    if new_bt is None:
        op_log.append(f"change_bond_order: unsupported order {new_order}")
        return False

    old_order = bond.GetBondTypeAsDouble()
    bond.SetBondType(new_bt)
    bond.SetIsAromatic(False)

    # Adjust H counts for valence
    order_diff = int(new_order) - int(old_order)
    if order_diff > 0:
        for idx in (i, j):
            atom = rw.GetAtomWithIdx(idx)
            atom.SetIsAromatic(False)
            atom.SetNoImplicit(False)
            total_hs = atom.GetTotalNumHs()
            reduce = min(order_diff, total_hs)
            if reduce > 0:
                atom.SetNumExplicitHs(max(0, total_hs - reduce))
                atom.SetNoImplicit(True)

    op_log.append(f"change_bond_order: {i}-{j} → order {new_order}")
    return True


def _op_ring_open(
    rw: RWMol,
    params: Dict[str, Any],
    op_log: List[str],
) -> bool:
    """Handle 'ring_open' operation: break a ring bond and add groups."""
    bond_pair = params.get("bond") or params.get("atoms")
    if not bond_pair or len(bond_pair) != 2:
        op_log.append("ring_open: missing 'bond' param")
        return False

    i, j = int(bond_pair[0]), int(bond_pair[1])
    bond = rw.GetBondBetweenAtoms(i, j)
    if bond is None:
        op_log.append(f"ring_open: no bond between {i} and {j}")
        return False

    rw.RemoveBond(i, j)

    # Add groups to the broken ends (default: add H)
    group_a = params.get("group_a", "H")
    group_b = params.get("group_b", "H")

    for idx, group in [(i, group_a), (j, group_b)]:
        if group == "H" or group is None:
            atom = rw.GetAtomWithIdx(idx)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            atom.SetNoImplicit(False)
        else:
            # Resolve and attach fragment
            frag_smi = FRAGMENT_TEMPLATES.get(group, group)
            frag_mol = Chem.MolFromSmiles(frag_smi)
            if frag_mol is None:
                op_log.append(f"ring_open: cannot parse group '{group}'")
                atom = rw.GetAtomWithIdx(idx)
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetNoImplicit(False)
                continue

            offset = rw.GetNumAtoms()
            combined = Chem.CombineMols(rw.GetMol(), frag_mol)
            rw_tmp = RWMol(combined)
            rw_tmp.AddBond(idx, offset, Chem.BondType.SINGLE)
            rw._replacement = rw_tmp.GetMol()  # type: ignore[attr-defined]
            # Re-wrap for subsequent operations
            rw = RWMol(rw._replacement)  # type: ignore[attr-defined]

    op_log.append(f"ring_open: broke ring bond {i}-{j}, groups=({group_a}, {group_b})")
    return True


def _op_modify_atom(
    rw: RWMol,
    params: Dict[str, Any],
    op_log: List[str],
) -> bool:
    """Handle 'modify_atom' operation: change element, charge, or H count."""
    atom_idx = params.get("atom_idx")
    if atom_idx is None:
        op_log.append("modify_atom: missing 'atom_idx'")
        return False

    atom_idx = int(atom_idx)
    if atom_idx >= rw.GetNumAtoms():
        op_log.append(f"modify_atom: atom_idx {atom_idx} out of range")
        return False

    atom = rw.GetAtomWithIdx(atom_idx)
    changes: List[str] = []

    if "element" in params:
        new_num = Chem.GetPeriodicTable().GetAtomicNumber(params["element"])
        atom.SetAtomicNum(new_num)
        changes.append(f"element→{params['element']}")

    if "charge" in params:
        atom.SetFormalCharge(int(params["charge"]))
        changes.append(f"charge→{params['charge']}")

    if "num_explicit_hs" in params:
        atom.SetNumExplicitHs(int(params["num_explicit_hs"]))
        changes.append(f"Hs→{params['num_explicit_hs']}")

    op_log.append(f"modify_atom: atom {atom_idx} [{', '.join(changes)}]")
    return True


def _op_fused_ring_retro(
    mol: Chem.Mol,
    params: Dict[str, Any],
    op_log: List[str],
) -> Optional[BreakResult]:
    """Handle 'fused_ring_retro' operation by delegating to the fused ring
    retro logic (FragmentOnBonds + capping).

    This is a terminal operation — it produces the final BreakResult directly.
    """
    bonds_to_break = params.get("bonds", [])
    if not bonds_to_break:
        op_log.append("fused_ring_retro: missing 'bonds'")
        return BreakResult(success=False, operation_log=op_log,
                           error="fused_ring_retro: no bonds specified")

    caps = params.get("caps", {})
    bond_order_changes = params.get("bond_order_changes", [])

    # Apply bond order changes first
    rw = RWMol(mol)
    for change in bond_order_changes:
        bp = change.get("bond", [])
        new_order = change.get("new_order", 1)
        if len(bp) != 2:
            continue
        ci, cj = int(bp[0]), int(bp[1])
        bond = rw.GetBondBetweenAtoms(ci, cj)
        if bond is None:
            continue
        order_map = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
                     3: Chem.BondType.TRIPLE}
        new_bt = order_map.get(int(new_order))
        if new_bt:
            old_order = bond.GetBondTypeAsDouble()
            bond.SetBondType(new_bt)
            bond.SetIsAromatic(False)
            diff = int(new_order) - int(old_order)
            if diff > 0:
                for aidx in (ci, cj):
                    a = rw.GetAtomWithIdx(aidx)
                    a.SetIsAromatic(False)
                    a.SetNoImplicit(False)
                    ths = a.GetTotalNumHs()
                    red = min(diff, ths)
                    if red > 0:
                        a.SetNumExplicitHs(max(0, ths - red))
                        a.SetNoImplicit(True)
            op_log.append(f"fused_ring_retro: bond_order {ci}-{cj} → {new_order}")

    try:
        Chem.SanitizeMol(rw)
    except Exception:
        try:
            Chem.SanitizeMol(
                rw,
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
            )
        except Exception:
            pass

    modified_mol = rw.GetMol()

    # Kekulize for aromatic systems
    try:
        Chem.Kekulize(modified_mol, clearAromaticFlags=True)
    except Exception:
        pass

    # Collect bond indices
    bond_indices = []
    for ai, aj in bonds_to_break:
        ai, aj = int(ai), int(aj)
        b = modified_mol.GetBondBetweenAtoms(ai, aj)
        if b is None:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"fused_ring_retro: bond ({ai}, {aj}) not found",
            )
        bond_indices.append(b.GetIdx())
        op_log.append(f"fused_ring_retro: marking bond {ai}-{aj}")

    # FragmentOnBonds
    dummy_labels = [(2 * k + 1, 2 * k + 2) for k in range(len(bond_indices))]
    try:
        fragmented = Chem.FragmentOnBonds(
            modified_mol, bond_indices, dummyLabels=dummy_labels,
        )
    except Exception as exc:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"FragmentOnBonds failed: {exc}",
        )

    if fragmented is None:
        return BreakResult(
            success=False, operation_log=op_log,
            error="FragmentOnBonds returned None",
        )

    # Split and process fragments
    frag_smiles_raw = Chem.MolToSmiles(fragmented)
    fragments = frag_smiles_raw.split(".")

    precursors: List[str] = []
    for frag_str in fragments:
        frag_str = frag_str.strip()
        if not frag_str:
            continue
        frag_mol = Chem.MolFromSmiles(frag_str)
        if frag_mol is None:
            continue

        # Replace dummy atoms: apply caps or default to H
        result_mol = frag_mol
        dummies = [
            (a.GetIdx(), a.GetIsotope())
            for a in result_mol.GetAtoms()
            if a.GetAtomicNum() == 0
        ]
        for _didx, isotope in sorted(dummies, key=lambda x: x[0], reverse=True):
            bond_k = (isotope - 1) // 2
            is_first = isotope % 2 == 1
            if bond_k < len(bonds_to_break):
                orig_atom = bonds_to_break[bond_k][0 if is_first else 1]
                side = "left" if is_first else "right"
            else:
                orig_atom = -1
                side = "left"

            # Look up cap
            cap = None
            for cap_key in [
                f"{orig_atom}:{side}",
                str(orig_atom),
                f"isotope_{isotope}",
                f"default_{side}",
                "default",
            ]:
                if cap_key in caps:
                    cap = caps[cap_key]
                    break

            # Apply cap
            result_mol = _apply_cap_to_dummy(result_mol, isotope, cap)
            if result_mol is None:
                break

        if result_mol is None:
            continue

        # Sanitize
        ok, _err = _sanitize_mol_fallback(result_mol)
        if not ok:
            continue
        try:
            smi = Chem.MolToSmiles(result_mol, canonical=True)
            if smi:
                precursors.append(smi)
        except Exception:
            continue

    if not precursors:
        return BreakResult(
            success=False, operation_log=op_log,
            error="fused_ring_retro: no valid precursors generated",
        )

    op_log.append(f"fused_ring_retro: produced {precursors}")
    return BreakResult(success=True, precursors=precursors, operation_log=op_log)


def _apply_cap_to_dummy(
    mol: Chem.Mol,
    isotope: int,
    cap: Optional[str],
) -> Optional[Chem.Mol]:
    """Replace a single dummy atom (by isotope) with *cap* or H."""
    dummy_idx = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == isotope:
            dummy_idx = atom.GetIdx()
            break
    if dummy_idx is None:
        return mol  # already processed

    if cap is None or cap == "H":
        # Remove dummy → effectively add H
        try:
            rw = RWMol(mol)
            rw.RemoveAtom(dummy_idx)
            return rw.GetMol()
        except Exception:
            return None

    if cap == "=O":
        # Special: replace dummy with double-bond oxygen
        try:
            rw = RWMol(mol)
            dummy_atom = rw.GetAtomWithIdx(dummy_idx)
            neighbors = [n.GetIdx() for n in dummy_atom.GetNeighbors()]
            if not neighbors:
                return None
            dummy_atom.SetAtomicNum(8)
            dummy_atom.SetIsotope(0)
            dummy_atom.SetFormalCharge(0)
            dummy_atom.SetNumExplicitHs(0)
            dummy_atom.SetNoImplicit(True)
            bond = rw.GetBondBetweenAtoms(dummy_idx, neighbors[0])
            if bond is not None:
                bond.SetBondType(Chem.BondType.DOUBLE)
            return rw.GetMol()
        except Exception:
            return None

    # General cap: use ReplaceSubstructs
    dummy_pat = Chem.MolFromSmarts("[#0]")
    replacement = Chem.MolFromSmiles(cap)
    if replacement is None:
        replacement = Chem.MolFromSmiles(f"[{cap}]")
    if replacement is None or dummy_pat is None:
        return None
    try:
        results = AllChem.ReplaceSubstructs(mol, dummy_pat, replacement, replaceAll=False)
        if results:
            return results[0]
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Public API: execute_operations
# ---------------------------------------------------------------------------

def execute_operations(
    smiles: str,
    operations: List[Dict[str, Any]],
    reaction_category: Optional[str] = None,
) -> BreakResult:
    """Execute a sequence of operations on a molecule.

    Each operation is a dict with a ``"type"`` key and type-specific
    parameters.  Supported types:

    - ``disconnect`` — Remove a bond (params: ``bond=[i,j]``)
    - ``add_fragment`` — Attach a fragment (params: ``atom_idx``,
      ``fragment``, optional ``attach_idx``, ``bond_type``)
    - ``change_bond_order`` — Change bond order (params: ``bond=[i,j]``,
      ``new_order``)
    - ``ring_open`` — Break a ring bond with optional groups
      (params: ``bond=[i,j]``, ``group_a``, ``group_b``)
    - ``modify_atom`` — Change atom properties (params: ``atom_idx``,
      optional ``element``, ``charge``, ``num_explicit_hs``)
    - ``fused_ring_retro`` — Fused ring retrosynthesis via
      ``FragmentOnBonds`` (params: ``bonds``, ``caps``,
      ``bond_order_changes``)

    Operations are executed sequentially on an ``RWMol``.  After all
    operations, the molecule is split via ``GetMolFrags`` and each
    fragment is sanitized.

    Parameters
    ----------
    smiles : str
        Target molecule SMILES.
    operations : List[Dict[str, Any]]
        Ordered list of operation dicts.
    reaction_category : Optional[str]
        Reaction category hint (currently unused, reserved for future).

    Returns
    -------
    BreakResult
        Contains ``success``, ``precursors``, ``operation_log``, ``error``.
    """
    op_log: List[str] = []

    # Validate input
    mol = parse_mol(smiles)
    if mol is None:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"invalid SMILES: {smiles}",
        )

    can = canonical(smiles)
    op_log.append(f"Target: {can}")
    if reaction_category:
        op_log.append(f"reaction_category: {reaction_category}")

    if not operations:
        return BreakResult(
            success=False, operation_log=op_log,
            error="no operations provided",
        )

    # Work on RWMol
    rw = RWMol(mol)

    for idx, op in enumerate(operations):
        op_type = op.get("type", "")
        op_log.append(f"Op {idx}: {op_type}")

        if op_type == "fused_ring_retro":
            # Terminal operation — produces BreakResult directly
            current_mol = rw.GetMol()
            result = _op_fused_ring_retro(current_mol, op, op_log)
            return result if result is not None else BreakResult(
                success=False, operation_log=op_log,
                error="fused_ring_retro failed",
            )

        elif op_type == "disconnect":
            ok = _op_disconnect(rw, op, op_log)

        elif op_type == "add_fragment":
            ok = _op_add_fragment(rw, op, op_log)
            # Check if replacement mol was created
            if ok and hasattr(rw, "_replacement"):
                rw = RWMol(rw._replacement)  # type: ignore[attr-defined]

        elif op_type == "change_bond_order":
            ok = _op_change_bond_order(rw, op, op_log)

        elif op_type == "ring_open":
            ok = _op_ring_open(rw, op, op_log)
            if ok and hasattr(rw, "_replacement"):
                rw = RWMol(rw._replacement)  # type: ignore[attr-defined]

        elif op_type == "modify_atom":
            ok = _op_modify_atom(rw, op, op_log)

        else:
            op_log.append(f"unknown operation type: {op_type}")
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"unknown operation type: {op_type}",
            )

        if not ok:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"operation {idx} ({op_type}) failed",
            )

    # --- Final: get mol, split fragments, sanitize ---
    try:
        edited_mol = rw.GetMol()
    except Exception as exc:
        return BreakResult(
            success=False, operation_log=op_log,
            error=f"failed to get mol after operations: {exc}",
        )

    # Three-level sanitize
    ok_san, err_san = _sanitize_mol_fallback(edited_mol)
    if not ok_san:
        return BreakResult(
            success=False, operation_log=op_log,
            error=err_san,
        )

    # Split into fragments
    try:
        frag_mols = Chem.GetMolFrags(edited_mol, asMols=True, sanitizeFrags=True)
    except Exception:
        try:
            frag_mols = Chem.GetMolFrags(edited_mol, asMols=True, sanitizeFrags=False)
        except Exception as exc:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"GetMolFrags failed: {exc}",
            )

    if not frag_mols:
        return BreakResult(
            success=False, operation_log=op_log,
            error="no fragments generated after operations",
        )

    precursors: List[str] = []
    for frag in frag_mols:
        ok_f, err_f = _sanitize_mol_fallback(frag)
        if not ok_f:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"fragment sanitize failed: {err_f}",
            )
        try:
            smi = Chem.MolToSmiles(frag, canonical=True)
        except Exception as exc:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"fragment to SMILES failed: {exc}",
            )
        if not smi:
            return BreakResult(
                success=False, operation_log=op_log,
                error="fragment produced empty SMILES",
            )
        # Verify round-trip
        check = Chem.MolFromSmiles(smi)
        if check is None:
            return BreakResult(
                success=False, operation_log=op_log,
                error=f"generated SMILES invalid: {smi}",
            )
        precursors.append(smi)

    op_log.append(f"execute_operations produced: {precursors}")
    return BreakResult(
        success=True,
        precursors=precursors,
        operation_log=op_log,
    )
