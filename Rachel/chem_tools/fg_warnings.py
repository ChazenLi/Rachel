"""M7: 官能团排斥与保护警告 — 检测官能团冲突、反应条件不兼容、保护基需求。"""

from __future__ import annotations

from typing import Any, Dict, List

from Rachel.display_language import dangerous_combo_text, selectivity_conflict_text
from Rachel.knowledge import get_base_profile

from ._rdkit_utils import load_template, parse_mol, smarts_match
from .fg_detect import detect_functional_groups, detect_protecting_groups

# ---------------------------------------------------------------------------
# Internal: reaction_category → reaction_family mapping
# ---------------------------------------------------------------------------
# Maps user-facing reaction category strings to the reaction family keys used
# in fg_compatibility_matrix.json.  Mirrors the mapping in forward_validate.py
# but kept local to avoid a cross-module dependency on M5.

_CATEGORY_FAMILY_MAP: Dict[str, List[str]] = {
    "oxidation": ["strong_oxidation", "mild_oxidation"],
    "strong_oxidation": ["strong_oxidation"],
    "mild_oxidation": ["mild_oxidation"],
    "reduction": ["strong_reduction", "mild_reduction"],
    "strong_reduction": ["strong_reduction"],
    "mild_reduction": ["mild_reduction"],
    "halogenation": ["radical"],
    "cross_coupling": ["pd_catalysis"],
    "coupling": ["pd_catalysis"],
    "amide_coupling": ["nucleophilic_addition"],
    "peptide_coupling": ["nucleophilic_addition"],
    "suzuki": ["pd_catalysis"],
    "heck": ["pd_catalysis"],
    "sonogashira": ["pd_catalysis"],
    "negishi": ["pd_catalysis"],
    "stille": ["pd_catalysis"],
    "buchwald": ["pd_catalysis"],
    "grignard": ["strong_base", "nucleophilic_addition"],
    "alkylation": ["strong_base"],
    "elimination": ["strong_base"],
    "acidic_ester_hydrolysis": ["strong_acid"],
    "base_promoted_ester_hydrolysis": ["strong_base"],
    "ester_hydrolysis": ["strong_acid", "strong_base"],
    "amide_formation": ["nucleophilic_addition"],
    "condensation": ["strong_acid"],
    "diazotization": ["diazotization"],
    "radical": ["radical"],
    "metathesis": ["metathesis"],
    "olefination": ["olefination"],
    "wittig": ["olefination"],
    "electrophilic_aromatic": ["electrophilic_aromatic"],
    "friedel_crafts": ["electrophilic_aromatic", "strong_acid"],
    "pd_catalysis": ["pd_catalysis"],
    "strong_base": ["strong_base"],
    "strong_acid": ["strong_acid"],
    "nucleophilic_addition": ["nucleophilic_addition"],
    "ir_catalysis": ["ir_catalysis"],
    "photoredox": ["photoredox_radical"],
    "photoredox_radical": ["photoredox_radical"],
    "tempo_oxidation": ["tempo_oxidation"],
    "organolithium": ["organolithium_cuprate"],
    "cuprate": ["organolithium_cuprate"],
    "organolithium_cuprate": ["organolithium_cuprate"],
    "fluoride_deprotection": ["fluoride_deprotection"],
}

_CONDITION_AMBIGUOUS_CATEGORY_TOKENS = {"oxidation", "reduction"}
_CONDITION_FALLBACK_CATEGORY_TOKENS = {"coupling"}


def _category_matches(category: str, key: str) -> bool:
    if category == key:
        return True
    category_tokens = [token for token in category.split("_") if token]
    key_tokens = [token for token in key.split("_") if token]
    if not category_tokens or not key_tokens:
        return False
    window = len(key_tokens)
    return any(
        category_tokens[index : index + window] == key_tokens
        for index in range(len(category_tokens) - window + 1)
    )


def resolve_reaction_condition_families(
    reaction_category: str,
    compatibility_matrix: Dict[str, Any] | None = None,
) -> List[str]:
    """Resolve only condition families justified by an explicit category."""
    category = reaction_category.lower().strip().replace("-", "_").replace(" ", "_")
    if category in _CATEGORY_FAMILY_MAP:
        families = list(_CATEGORY_FAMILY_MAP[category])
    else:
        matched: List[tuple[str, List[str]]] = []
        for key, mapped in _CATEGORY_FAMILY_MAP.items():
            if key in _CONDITION_AMBIGUOUS_CATEGORY_TOKENS:
                continue
            if _category_matches(category, key):
                matched.append((key, mapped))
        if any(key not in _CONDITION_FALLBACK_CATEGORY_TOKENS for key, _ in matched):
            matched = [
                (key, mapped)
                for key, mapped in matched
                if key not in _CONDITION_FALLBACK_CATEGORY_TOKENS
            ]
        families = [family for _, mapped in matched for family in mapped]
        if not families and compatibility_matrix and category in compatibility_matrix:
            families = [category]

    return list(dict.fromkeys(families))


# ---------------------------------------------------------------------------
# Internal: superclass grouping for competing-groups detection
# ---------------------------------------------------------------------------

_SUPERCLASS: Dict[str, str] = {
    "primary_alcohol": "alcohol",
    "secondary_alcohol": "alcohol",
    "tertiary_alcohol": "alcohol",
    "benzylic_alcohol": "alcohol",
    "allylic_alcohol": "alcohol",
    "phenol": "phenol",
    "primary_amine": "amine",
    "secondary_amine": "amine",
    "aromatic_amine": "amine",
    "aldehyde": "carbonyl",
    "ketone": "carbonyl",
    "cyclic_ketone": "carbonyl",
    "carboxylic_acid": "acid",
    "alkene": "unsaturated",
    "terminal_alkene": "unsaturated",
    "alkyne": "unsaturated",
    "terminal_alkyne": "unsaturated",
    "ester_generic": "ester",
    "lactone": "ester",
    "thiol": "thiol",
    "aryl_halide": "halide",
    "aryl_bromide": "halide",
    "aryl_iodide": "halide",
    "aryl_chloride": "halide",
}


# ---------------------------------------------------------------------------
# check_fg_conflicts
# ---------------------------------------------------------------------------

def check_fg_conflicts(smiles: str, knowledge_profile=None) -> Dict[str, Any]:
    """Detect functional-group conflicts within a molecule.

    Calls :func:`detect_functional_groups` to obtain the FGs present, then
    loads ``selectivity_conflicts.json`` and ``dangerous_combos.json`` to
    identify pairwise conflicts, competing groups, and dangerous combinations.

    Returns
    -------
    dict
        ``{"ok": True, "conflicts": [...], "competing_groups": [...],
          "dangerous_combos": [...]}``
        or ``{"ok": False, "error": "...", "input": smiles}`` on invalid input.
    """
    # --- validate input ---
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    fg_result = detect_functional_groups(smiles, knowledge_profile=knowledge_profile)
    if not fg_result.get("ok"):
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    detected_groups = fg_result.get("groups", {})
    detected_names = set(detected_groups.keys())

    # --- load templates ---
    conflicts_data = load_template(
        "selectivity_conflicts.json", knowledge_profile=knowledge_profile
    )
    dangerous_data = load_template(
        "dangerous_combos.json", knowledge_profile=knowledge_profile
    )
    reactivity_data = load_template(
        "selectivity_reactivity.json", knowledge_profile=knowledge_profile
    )

    # Build a reverse map: fg_key → set of selectivity categories
    fg_key_to_categories: Dict[str, set] = {}
    for cat, info in reactivity_data.items():
        for fk in info.get("fg_keys", []):
            fg_key_to_categories.setdefault(fk, set()).add(cat)

    # Determine which selectivity categories are present
    present_categories: set = set()
    for fg_name in detected_names:
        for cat_name in fg_key_to_categories.get(fg_name, set()):
            present_categories.add(cat_name)

    # --- 1. Pairwise conflicts from selectivity_conflicts.json ---
    conflict_list: List[Dict[str, Any]] = []
    pairs = conflicts_data.get("pairs", [])
    for entry in pairs:
        if len(entry) < 5:
            continue
        cat_a, cat_b, reaction_type, severity, note = (
            entry[0], entry[1], entry[2], entry[3], entry[4]
        )
        if cat_a in present_categories and cat_b in present_categories:
            conflict_list.append({
                "group_a": cat_a,
                "group_b": cat_b,
                "reaction_type": reaction_type,
                "severity": severity,
                "note": selectivity_conflict_text(
                    cat_a, cat_b, reaction_type, note
                ),
            })

    # --- 2. Competing groups (same superclass present multiple times) ---
    competing_groups: List[str] = []
    superclass_members: Dict[str, List[str]] = {}
    for fg_name in detected_names:
        sc = _SUPERCLASS.get(fg_name)
        if sc:
            superclass_members.setdefault(sc, []).append(fg_name)

    for sc, members in superclass_members.items():
        if len(members) >= 2:
            competing_groups.extend(members)

    # Also detect same FG appearing multiple times (count >= 2)
    for fg_name, info in detected_groups.items():
        if info.get("count", 0) >= 2 and fg_name not in competing_groups:
            competing_groups.append(fg_name)

    # --- 3. Dangerous combos from dangerous_combos.json ---
    dangerous_list: List[Dict[str, Any]] = []
    for combo_name, combo_data in dangerous_data.items():
        if not isinstance(combo_data, list) or len(combo_data) < 3:
            continue
        smarts_a, smarts_b, warning = combo_data[0], combo_data[1], combo_data[2]
        match_a = smarts_match(mol, smarts_a)
        match_b = smarts_match(mol, smarts_b)
        if match_a and match_b:
            dangerous_list.append({
                "groups": [combo_name],
                "warning": dangerous_combo_text(combo_name, warning),
            })

    return {
        "ok": len(conflict_list) == 0 and len(dangerous_list) == 0,
        "conflicts": conflict_list,
        "competing_groups": competing_groups,
        "dangerous_combos": dangerous_list,
    }


# ---------------------------------------------------------------------------
# check_reaction_compatibility
# ---------------------------------------------------------------------------

def check_reaction_compatibility(
    smiles: str, reaction_category: str, knowledge_profile=None
) -> Dict[str, Any]:
    """Check functional-group compatibility with a given reaction category.

    Loads ``fg_compatibility_matrix.json`` and checks every detected FG
    against the reaction families associated with *reaction_category*.

    Returns
    -------
    dict
        ``{"compatible": True/False, "warnings": [...]}``
        where each warning is ``{"fg": str, "level": str,
        "reaction_family": str, "note": str}``.
        ``level`` is one of ``"risky"`` or ``"forbidden"``.
        ``compatible`` is ``False`` when any ``"forbidden"`` warning exists.

        On invalid SMILES: ``{"ok": False, "error": "...", "input": smiles}``.
    """
    # --- validate input ---
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    profile = knowledge_profile or get_base_profile()

    fg_result = detect_functional_groups(smiles, knowledge_profile=profile)
    if not fg_result.get("ok"):
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    detected_groups = fg_result.get("groups", {})
    detected_names = set(detected_groups.keys())

    # --- load compatibility matrix ---
    compat_matrix = load_template(
        "fg_compatibility_matrix.json", knowledge_profile=profile
    )

    # --- resolve reaction_category to reaction families ---
    families = resolve_reaction_condition_families(reaction_category, compat_matrix)

    # --- check each detected FG against each reaction family ---
    warnings: List[Dict[str, Any]] = []
    has_forbidden = False
    matched_families: set[str] = set()

    for family in families:
        family_rules = compat_matrix.get(family, {})
        applies_to = family_rules.get("__applies_to", "")

        for fg_name in detected_names:
            level = family_rules.get(fg_name)
            if level is None or level not in ("risky", "forbidden"):
                continue

            if level == "forbidden":
                has_forbidden = True
            matched_families.add(family)

            note = (
                f"{fg_name} is {level} under {family} conditions"
                f" ({applies_to})" if applies_to else
                f"{fg_name} is {level} under {family} conditions"
            )

            warnings.append({
                "fg": fg_name,
                "level": level,
                "reaction_family": family,
                "note": note,
            })

    return {
        "compatible": not has_forbidden,
        "warnings": warnings,
        "knowledge_profile_hash": profile.digest,
        "knowledge_refs": [
            profile.source("chem.fg_compatibility", family)
            for family in sorted(matched_families)
        ],
    }


# ---------------------------------------------------------------------------
# Internal: FG name → protectable type mapping
# ---------------------------------------------------------------------------
# Maps functional group names (from functional_groups.json) to the protection
# type categories used in protecting_groups.json definitions.

_FG_TO_PROTECTION_TYPE: Dict[str, str] = {
    # amine-type FGs
    "primary_amine": "amine",
    "secondary_amine": "amine",
    "aromatic_amine": "amine",
    # hydroxyl-type FGs
    "primary_alcohol": "hydroxyl",
    "secondary_alcohol": "hydroxyl",
    "tertiary_alcohol": "hydroxyl",
    "benzylic_alcohol": "hydroxyl",
    "allylic_alcohol": "hydroxyl",
    "phenol": "hydroxyl",
    # carboxyl-type FGs
    "carboxylic_acid": "carboxyl",
    # thiol-type FGs
    "thiol": "thiol",
    # carbonyl-type FGs
    "aldehyde": "carbonyl",
    "ketone": "carbonyl",
    "cyclic_ketone": "carbonyl",
}

# ---------------------------------------------------------------------------
# Internal: default recommended PGs per protection type (ordered by preference)
# ---------------------------------------------------------------------------

_DEFAULT_PG_BY_TYPE: Dict[str, List[str]] = {
    "amine": ["Boc", "Fmoc", "Cbz"],
    "hydroxyl": ["TBS", "Bn_O", "Acetyl_O"],
    "carboxyl": ["tBuEster", "MethylEster", "BnEster"],
    "thiol": ["Trityl_S", "Acm_S", "StBu_S"],
    "carbonyl": ["Acetal", "Dioxolane", "Dithiane"],
}

# ---------------------------------------------------------------------------
# Internal: removal condition keywords → reaction families
# ---------------------------------------------------------------------------
# Used by check_deprotection_safety to infer which reaction families are
# involved in a deprotection step, based on the "removal" text in
# protecting_groups.json.

_REMOVAL_KEYWORD_TO_FAMILIES: List[tuple] = [
    # (keyword_in_removal_text, [reaction_families])
    ("TFA", ["strong_acid"]),
    ("HCl", ["strong_acid"]),
    ("HF", ["strong_acid"]),
    ("TCA", ["strong_acid"]),
    ("acid-labile", ["strong_acid"]),
    ("dilute acid", ["strong_acid"]),
    ("PPTS", ["strong_acid"]),
    ("p-TsOH", ["strong_acid"]),
    ("AcOH", ["strong_acid"]),
    ("NaOH", ["strong_base"]),
    ("LiOH", ["strong_base"]),
    ("K2CO3", ["strong_base"]),
    ("piperidine", ["strong_base"]),
    ("DBU", ["strong_base"]),
    ("Et3N", ["strong_base"]),
    ("base-labile", ["strong_base"]),
    ("mild base", ["strong_base"]),
    ("TBAF", ["fluoride_deprotection"]),
    ("fluoride", ["fluoride_deprotection"]),
    ("H2/Pd", ["mild_reduction"]),
    ("hydrogenolysis", ["mild_reduction"]),
    ("LiAlH4", ["strong_reduction"]),
    ("DIBAL", ["strong_reduction"]),
    ("Na/naphthalene", ["strong_reduction"]),
    ("SmI2", ["strong_reduction"]),
    ("Zn/AcOH", ["strong_reduction"]),
    ("reductive", ["strong_reduction"]),
    ("DDQ", ["mild_oxidation"]),
    ("CAN", ["mild_oxidation"]),
    ("oxidative", ["mild_oxidation"]),
    ("NBS", ["radical"]),
    ("I2", ["mild_oxidation"]),
    ("Hg(OAc)2", ["strong_acid"]),
    ("Pd(0)", ["pd_catalysis"]),
    ("N2H4", ["nucleophilic_addition"]),
    ("hydrazinolysis", ["nucleophilic_addition"]),
    ("PhSH", ["nucleophilic_addition"]),
    ("thiolysis", ["nucleophilic_addition"]),
    ("thiourea", ["nucleophilic_addition"]),
]


# ---------------------------------------------------------------------------
# suggest_protection_needs
# ---------------------------------------------------------------------------

def suggest_protection_needs(
    smiles: str, planned_reaction: str, knowledge_profile=None
) -> Dict[str, Any]:
    """Suggest functional groups that need protection for a planned reaction.

    Calls :func:`detect_functional_groups` and :func:`detect_protecting_groups`
    to obtain present FGs and existing PGs, then uses
    ``fg_compatibility_matrix.json`` to find incompatible FGs under the
    planned reaction conditions.  For each incompatible FG, recommends a
    suitable protecting group from ``protecting_groups.json`` and checks
    orthogonality with existing PGs.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.
    planned_reaction : str
        Reaction category string (e.g. ``"grignard"``, ``"oxidation"``).

    Returns
    -------
    dict
        ``{"needs_protection": bool, "suggestions": [...],
          "existing_pg": [...], "orthogonality_notes": [...]}``
        or ``{"ok": False, "error": "...", "input": smiles}`` on invalid input.
    """
    # --- validate input ---
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    fg_result = detect_functional_groups(smiles, knowledge_profile=knowledge_profile)
    if not fg_result.get("ok"):
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    pg_result = detect_protecting_groups(smiles, knowledge_profile=knowledge_profile)

    detected_groups = fg_result.get("groups", {})
    detected_names = set(detected_groups.keys())

    # Collect existing PG names
    existing_pg: List[str] = []
    if pg_result.get("status") == "ok":
        for pg_entry in pg_result.get("detected", []):
            existing_pg.append(pg_entry.get("name", ""))

    # --- resolve planned_reaction to reaction families ---
    # --- load templates ---
    compat_matrix = load_template(
        "fg_compatibility_matrix.json", knowledge_profile=knowledge_profile
    )
    families = resolve_reaction_condition_families(planned_reaction, compat_matrix)
    pg_data = load_template(
        "protecting_groups.json", knowledge_profile=knowledge_profile
    )
    pg_defs = {
        k: v for k, v in pg_data.get("definitions", {}).items()
        if not k.startswith("__")
    }
    orthogonal_map = pg_data.get("orthogonal_map", {})

    # --- find incompatible FGs ---
    suggestions: List[Dict[str, Any]] = []
    for family in families:
        family_rules = compat_matrix.get(family, {})
        for fg_name in detected_names:
            level = family_rules.get(fg_name)
            if level not in ("risky", "forbidden"):
                continue

            # Determine protection type for this FG
            prot_type = _FG_TO_PROTECTION_TYPE.get(fg_name)
            if prot_type is None:
                continue  # FG not protectable

            # Find a recommended PG
            candidates = _DEFAULT_PG_BY_TYPE.get(prot_type, [])
            recommended_pg = candidates[0] if candidates else None
            if recommended_pg is None:
                continue

            reason = (
                f"{fg_name} is {level} under {family} conditions; "
                f"recommend {recommended_pg} ({prot_type} protection)"
            )

            # Avoid duplicate suggestions for the same FG
            already = any(s["fg"] == fg_name for s in suggestions)
            if already:
                continue

            suggestions.append({
                "fg": fg_name,
                "recommended_pg": recommended_pg,
                "reason": reason,
            })

    # --- orthogonality check ---
    orthogonality_notes: List[str] = []
    suggested_pg_names = [s["recommended_pg"] for s in suggestions]
    all_pg_names = existing_pg + suggested_pg_names

    if len(all_pg_names) >= 2:
        checked: set = set()
        for i, pg_a in enumerate(all_pg_names):
            for pg_b in all_pg_names[i + 1:]:
                pair = tuple(sorted([pg_a, pg_b]))
                if pair in checked:
                    continue
                checked.add(pair)
                orth_a = set(orthogonal_map.get(pg_a, []))
                orth_b = set(orthogonal_map.get(pg_b, []))
                if pg_b in orth_a or pg_a in orth_b:
                    orthogonality_notes.append(
                        f"{pg_a} and {pg_b} are orthogonal — can be removed independently"
                    )
                else:
                    orthogonality_notes.append(
                        f"{pg_a} and {pg_b} may not be orthogonal — verify removal compatibility"
                    )

    return {
        "needs_protection": len(suggestions) > 0,
        "suggestions": suggestions,
        "existing_pg": existing_pg,
        "orthogonality_notes": orthogonality_notes,
    }


# ---------------------------------------------------------------------------
# check_deprotection_safety
# ---------------------------------------------------------------------------

def check_deprotection_safety(
    smiles: str, pg_to_remove: str, knowledge_profile=None
) -> Dict[str, Any]:
    """Check whether removing a protecting group is safe for other FGs.

    Looks up the deprotection conditions for *pg_to_remove* in
    ``protecting_groups.json``, infers the reaction families involved, then
    checks every detected functional group against those families using
    ``fg_compatibility_matrix.json``.

    Parameters
    ----------
    smiles : str
        Molecule SMILES.
    pg_to_remove : str
        Name of the protecting group to remove (e.g. ``"Boc"``, ``"TBS"``).

    Returns
    -------
    dict
        ``{"safe": bool, "risks": [...]}``
        where each risk is ``{"affected_fg": str, "condition_family": str,
        "level": str}``.

        On invalid SMILES: ``{"ok": False, "error": "...", "input": smiles}``.
    """
    # --- validate input ---
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    fg_result = detect_functional_groups(smiles, knowledge_profile=knowledge_profile)
    if not fg_result.get("ok"):
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    detected_groups = fg_result.get("groups", {})
    detected_names = set(detected_groups.keys())

    # --- load protecting_groups.json and find the PG entry ---
    pg_data = load_template(
        "protecting_groups.json", knowledge_profile=knowledge_profile
    )
    pg_defs = {
        k: v for k, v in pg_data.get("definitions", {}).items()
        if not k.startswith("__")
    }

    pg_info = pg_defs.get(pg_to_remove)
    if pg_info is None:
        return {
            "safe": True,
            "risks": [],
            "note": f"protecting group '{pg_to_remove}' not found in database",
        }

    removal_text = pg_info.get("removal", "")

    # --- infer reaction families from removal conditions ---
    deprotection_families: List[str] = []
    seen_fam: set = set()
    for keyword, fam_list in _REMOVAL_KEYWORD_TO_FAMILIES:
        if keyword.lower() in removal_text.lower():
            for fam in fam_list:
                if fam not in seen_fam:
                    seen_fam.add(fam)
                    deprotection_families.append(fam)

    # --- check FG compatibility against deprotection families ---
    compat_matrix = load_template(
        "fg_compatibility_matrix.json", knowledge_profile=knowledge_profile
    )

    risks: List[Dict[str, Any]] = []
    for family in deprotection_families:
        family_rules = compat_matrix.get(family, {})
        for fg_name in detected_names:
            level = family_rules.get(fg_name)
            if level not in ("risky", "forbidden"):
                continue
            risks.append({
                "affected_fg": fg_name,
                "condition_family": family,
                "level": level,
            })

    return {
        "safe": len(risks) == 0,
        "risks": risks,
    }
