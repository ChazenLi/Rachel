"""M2: 官能团识别与匹配 — 识别官能团、反应活性位点、保护基，提供官能团→反应映射。"""

from __future__ import annotations

from typing import Any, Dict

from ._rdkit_utils import canonical, load_template, parse_mol, smarts_match


# ---------------------------------------------------------------------------
# detect_functional_groups
# ---------------------------------------------------------------------------

def detect_functional_groups(smiles: str) -> Dict[str, Any]:
    """Identify all functional groups present in a molecule.

    Loads ``templates/functional_groups.json`` and runs each SMARTS pattern
    against the molecule.  For ``bridgehead_atom`` entries the matches are
    deduplicated by the centre atom index (first element of each match tuple)
    so that the same bridgehead is not counted multiple times from different
    directions.  Groups with zero hits are **not** included in the output.

    Returns
    -------
    dict
        ``{"ok": True, "smiles": <canonical>, "groups": {name: {"count": int, "atoms": [[int]]}}}``
        or ``{"ok": False, "error": "...", "input": smiles}`` on invalid input.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    can = canonical(smiles)
    fg_template = load_template("functional_groups.json")

    groups: Dict[str, Any] = {}
    for name, sma in fg_template.items():
        # Skip comment keys and non-string values
        if name.startswith("__comment") or not isinstance(sma, str):
            continue
        matches = smarts_match(mol, sma)
        if not matches:
            continue

        # bridgehead_atom: deduplicate by centre atom (first atom in tuple)
        if name == "bridgehead_atom":
            seen_centers: set[int] = set()
            unique: list = []
            for m in matches:
                center = m[0]
                if center not in seen_centers:
                    seen_centers.add(center)
                    unique.append(m)
            matches = tuple(unique)

        groups[name] = {
            "count": len(matches),
            "atoms": [list(m) for m in matches],
        }

    return {"ok": True, "smiles": can, "groups": groups}


# ---------------------------------------------------------------------------
# detect_reactive_sites
# ---------------------------------------------------------------------------

def detect_reactive_sites(smiles: str) -> Dict[str, Any]:
    """Identify reactive sites and group them by category prefix.

    Loads ``templates/reactive_sites.json``, matches each SMARTS, then groups
    the results by the prefix before the first underscore (e.g.
    ``"electrophilic_carbonyl"`` → category ``"electrophilic"``).

    Returns
    -------
    dict
        ``{"ok": True, "sites": {category: {site_name: {"count": int, "atoms": [[int]]}}}}``
        or ``{"ok": False, "error": "...", "input": smiles}`` on invalid input.
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES", "input": smiles}

    rs_template = load_template("reactive_sites.json")

    sites: Dict[str, Dict[str, Any]] = {}
    for name, sma in rs_template.items():
        # Skip comment keys and non-string values
        if name.startswith("__comment") or not isinstance(sma, str):
            continue
        matches = smarts_match(mol, sma)
        if not matches:
            continue

        parts = name.split("_")
        category = parts[0] if parts else "other"

        if category not in sites:
            sites[category] = {}

        sites[category][name] = {
            "count": len(matches),
            "atoms": [list(m) for m in matches],
        }

    return {"ok": True, "sites": sites}


# ---------------------------------------------------------------------------
# detect_protecting_groups
# ---------------------------------------------------------------------------

# Substructure containment: when the specific PG is detected, remove the
# generic one (e.g. TBS hit removes TMS because TBS contains TMS substructure).
_SUBSUMES: Dict[str, str] = {
    # Silyl ethers: TBS/TBDPS/TIPS/TES subsume TMS
    "TBS": "TMS",
    "TBDPS": "TMS",
    "TIPS": "TMS",
    "TES": "TMS",
    # Benzyl derivatives: PMB/DMB/Nap subsume Bn
    "PMB_O": "Bn_O",
    "DMB_O": "Bn_O",
    "Nap_O": "Bn_O",
    "PMB_N": "Bn_N",
    "Mob_S": "Bn_S",
    # Trityl derivatives
    "MMTr_O": "Trityl_O",
    "Mtt_N": "Trityl_O",
    "Mtt_S": "Trityl_S",
    # Chain ethers: MEM/SEM subsume MOM
    "MEM": "MOM",
    "SEM": "MOM",
    # Acyl: Lev subsumes Acetyl_O
    "Lev_O": "Acetyl_O",
    # Ketal/acetal: Acetonide subsumes Ketal, Dioxolane subsumes Acetal/Ketal,
    # Benzylidene subsumes Acetal
    "Acetonide": "Ketal",
    "Dioxolane": "Acetal",
    "Benzylidene": "Acetal",
}

_ALL_PG_TYPES = ("amine", "hydroxyl", "carboxyl", "thiol", "carbonyl", "phosphate")


def _deduplicate_pgs(
    detected: list,
) -> tuple:
    """Remove redundant PG detections based on substructure containment.

    Returns ``(deduped_list, log_messages)``.
    """
    from typing import Set, List as TList

    names_found: Set[str] = {d["name"] for d in detected}
    to_remove: Set[str] = set()
    logs: TList[str] = []

    for specific, generic in _SUBSUMES.items():
        if specific in names_found and generic in names_found:
            to_remove.add(generic)
            logs.append(f"dedup: {generic} subsumed by {specific}")

    if not to_remove:
        return detected, logs

    return [d for d in detected if d["name"] not in to_remove], logs


def _check_orthogonality(detected: list, orthogonal_map: dict) -> list:
    """Bidirectional orthogonality check for detected protecting groups."""
    notes: list = []
    if len(detected) < 2:
        return notes

    for i, d1 in enumerate(detected):
        for d2 in detected[i + 1 :]:
            n1, n2 = d1["name"], d2["name"]
            orth_1 = set(orthogonal_map.get(n1, []))
            orth_2 = set(orthogonal_map.get(n2, []))
            if n2 in orth_1 or n1 in orth_2:
                notes.append(f"{n1} and {n2} are orthogonal — can be removed independently")
            else:
                notes.append(
                    f"{n1} and {n2} may not be orthogonal — verify removal condition compatibility"
                )

    if len(detected) >= 3:
        notes.append(
            "Multiple protecting groups detected — plan removal order to avoid cross-reactivity"
        )

    return notes


def detect_protecting_groups(smiles: str) -> Dict[str, Any]:
    """Detect existing protecting groups in a molecule.

    Loads ``templates/protecting_groups.json``, precompiles SMARTS patterns,
    performs substructure containment deduplication (e.g. TBS hit removes TMS),
    and runs bidirectional orthogonality checks.

    Returns
    -------
    dict
        On success::

            {
                "status": "ok",
                "detected": [{"name", "type", "count", "removal", "orthogonal_to"}, ...],
                "summary": {"has_protected_amine": bool, ..., "total_pg_count": int},
                "orthogonality_notes": [str, ...]
            }

        On invalid SMILES::

            {"status": "invalid_smiles", "error": "...", "input": smiles}
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"status": "invalid_smiles", "error": "invalid SMILES", "input": smiles}

    pg_data = load_template("protecting_groups.json")
    pg_defs = pg_data.get("definitions", {})
    orthogonal_map = pg_data.get("orthogonal_map", {})

    # Precompile SMARTS and match
    from rdkit import Chem as _Chem

    detected: list = []
    for name, info in pg_defs.items():
        if not isinstance(info, dict) or "smarts" not in info:
            continue
        pat = _Chem.MolFromSmarts(info["smarts"])
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat, uniquify=True)
        if not matches:
            continue
        detected.append(
            {
                "name": name,
                "type": info.get("type", ""),
                "count": len(matches),
                "removal": info.get("removal", ""),
                "orthogonal_to": orthogonal_map.get(name, []),
            }
        )

    if not detected:
        return {
            "status": "no_pg_found",
            "detected": [],
            "summary": {f"has_protected_{t}": False for t in _ALL_PG_TYPES}
            | {"total_pg_count": 0},
            "orthogonality_notes": [],
        }

    # Substructure containment dedup
    detected, _dedup_logs = _deduplicate_pgs(detected)

    # Summary
    types_found = {d["type"] for d in detected}
    total = sum(d["count"] for d in detected)
    summary: Dict[str, Any] = {
        f"has_protected_{t}": t in types_found for t in _ALL_PG_TYPES
    }
    summary["total_pg_count"] = total

    # Orthogonality
    ortho_notes = _check_orthogonality(detected, orthogonal_map)

    return {
        "status": "ok",
        "detected": detected,
        "summary": summary,
        "orthogonality_notes": ortho_notes,
    }


# ---------------------------------------------------------------------------
# get_fg_reaction_mapping
# ---------------------------------------------------------------------------

# Built-in FG name → possible reaction types mapping table.
# Distilled from the old system's _TRANSFORMATION_LABELS and domain knowledge.
_FG_REACTION_MAP: Dict[str, list] = {
    # Carbonyl electrophiles
    "aldehyde": ["nucleophilic_addition", "wittig_olefination", "aldol_condensation", "reductive_amination", "oxidation", "reduction"],
    "ketone": ["nucleophilic_addition", "wittig_olefination", "reductive_amination", "baeyer_villiger", "reduction"],
    "cyclic_ketone": ["nucleophilic_addition", "wittig_olefination", "baeyer_villiger", "reduction", "ring_opening"],
    "alpha_diketone": ["nucleophilic_addition", "reduction", "oxidative_cleavage"],
    # Carboxylic acid derivatives
    "carboxylic_acid": ["esterification", "amide_coupling", "reduction", "decarboxylation"],
    "ester_generic": ["hydrolysis", "reduction", "transesterification", "claisen_condensation"],
    "lactone": ["hydrolysis", "reduction", "ring_opening"],
    "amide_generic": ["hydrolysis", "reduction"],
    "primary_amide": ["dehydration_to_nitrile", "hydrolysis", "reduction"],
    "secondary_amide": ["hydrolysis", "reduction"],
    "tertiary_amide": ["hydrolysis", "reduction"],
    "lactam": ["hydrolysis", "reduction", "ring_opening"],
    "beta_lactam": ["ring_opening", "hydrolysis"],
    "anhydride": ["acylation", "hydrolysis", "aminolysis"],
    "acyl_halide": ["acylation", "aminolysis", "reduction"],
    "acyl_chloride": ["acylation", "aminolysis", "reduction"],
    # Unsaturated electrophiles
    "alpha_beta_unsat_carbonyl": ["conjugate_addition", "reduction"],
    "alpha_beta_unsat_ester": ["conjugate_addition", "reduction"],
    "michael_acceptor": ["conjugate_addition"],
    # C=N electrophiles
    "imine": ["nucleophilic_addition", "reduction"],
    "nitrile": ["nucleophilic_addition", "hydrolysis", "reduction"],
    "isocyanate": ["addition_with_nucleophile"],
    # Epoxide / aziridine
    "epoxide": ["ring_opening", "reduction"],
    "aziridine": ["ring_opening"],
    # Halides
    "alkyl_halide": ["SN2_substitution", "elimination", "grignard_formation"],
    "alkyl_chloride": ["SN2_substitution", "elimination"],
    "alkyl_bromide": ["SN2_substitution", "elimination", "grignard_formation"],
    "alkyl_iodide": ["SN2_substitution", "elimination"],
    "benzyl_halide": ["SN2_substitution"],
    "allyl_halide": ["SN2_substitution", "allylic_substitution"],
    "aryl_halide": ["cross_coupling", "nucleophilic_aromatic_substitution"],
    "aryl_fluoride": ["nucleophilic_aromatic_substitution"],
    "aryl_chloride": ["cross_coupling", "buchwald_coupling"],
    "aryl_bromide": ["cross_coupling", "suzuki_coupling", "heck_coupling"],
    "aryl_iodide": ["cross_coupling", "suzuki_coupling", "heck_coupling", "sonogashira_coupling"],
    "vinyl_halide": ["cross_coupling", "elimination"],
    # Leaving groups
    "mesylate": ["SN2_substitution", "elimination"],
    "tosylate": ["SN2_substitution", "elimination"],
    "triflate": ["cross_coupling", "SN2_substitution"],
    # Alcohols / phenols
    "alcohol": ["oxidation", "protection", "esterification", "elimination"],
    "primary_alcohol": ["oxidation", "protection", "esterification"],
    "secondary_alcohol": ["oxidation", "protection", "esterification"],
    "tertiary_alcohol": ["elimination", "protection"],
    "phenol": ["protection", "ether_formation", "electrophilic_aromatic_substitution"],
    "diol": ["protection", "oxidative_cleavage", "acetonide_formation"],
    # Amines
    "primary_amine": ["protection", "reductive_amination", "amide_coupling", "diazotization"],
    "secondary_amine": ["protection", "reductive_amination", "amide_coupling"],
    "tertiary_amine": ["quaternization", "cope_elimination"],
    "aniline": ["protection", "diazotization", "electrophilic_aromatic_substitution"],
    # Alkenes / alkynes
    "alkene": ["hydrogenation", "epoxidation", "dihydroxylation", "hydroboration", "metathesis"],
    "terminal_alkene": ["hydrogenation", "hydroboration", "cross_metathesis", "wacker_oxidation"],
    "alkyne": ["hydrogenation", "hydration", "sonogashira_coupling"],
    "terminal_alkyne": ["sonogashira_coupling", "click_chemistry", "deprotonation"],
    # Sulfur
    "thiol": ["protection", "disulfide_formation", "thioether_formation"],
    "sulfide": ["oxidation_to_sulfoxide"],
    "sulfoxide": ["oxidation_to_sulfone", "pummerer_rearrangement"],
    "sulfone": ["julia_olefination", "elimination"],
    "sulfonyl_chloride": ["sulfonamide_formation"],
    # Nitro
    "nitro": ["reduction_to_amine"],
    "nitro_aromatic": ["reduction_to_amine", "nucleophilic_aromatic_substitution"],
    # Boronic acids
    "boronic_acid": ["suzuki_coupling"],
    "boronate_ester": ["suzuki_coupling"],
    # Phosphorus
    "phosphonate": ["horner_wadsworth_emmons"],
    # Acetal / ketal
    "acetal": ["hydrolysis_to_aldehyde"],
    "ketal": ["hydrolysis_to_ketone"],
    # Oxime / hydrazone
    "oxime": ["beckmann_rearrangement", "hydrolysis"],
    "hydrazone": ["wolff_kishner_reduction", "hydrolysis"],
    # Carbamate
    "carbamate": ["deprotection", "hydrolysis"],
}


def get_fg_reaction_mapping(groups: Dict) -> Dict[str, list]:
    """Return FG → possible reaction types for groups actually present.

    Parameters
    ----------
    groups : dict
        The ``"groups"`` value from :func:`detect_functional_groups` output,
        i.e. ``{fg_name: {"count": int, "atoms": [[int]]}}``.

    Returns
    -------
    Dict[str, List[str]]
        Only entries for FGs that exist in *groups* **and** have a known
        mapping are returned.
    """
    if not groups or not isinstance(groups, dict):
        return {}

    result: Dict[str, list] = {}
    for fg_name in groups:
        reactions = _FG_REACTION_MAP.get(fg_name)
        if reactions:
            result[fg_name] = list(reactions)  # copy to avoid mutation
    return result
