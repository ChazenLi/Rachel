"""M5: forward validation and atom balance."""

from __future__ import annotations

import logging
from collections import Counter
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions

from ._rdkit_utils import (
    canonical,
    load_template,
    mol_formula_counter,
    parse_mol,
    smarts_match,
    tanimoto,
)
from .atom_mapping_audit import audit_atom_mapping
from .eas_site_audit import audit_eas_site_selectivity
from .fg_delta_audit import audit_fg_delta
from .fg_detect import detect_functional_groups
from .fg_warnings import check_reaction_compatibility
from .graph_delta_audit import audit_graph_delta
from .reaction_family_validate import validate_reaction_family
from .ring_topology_audit import audit_ring_topology
from .site_audit import audit_site_retention
from .validation_evidence_packet import build_validation_evidence_packet
from .validation_policy import build_validation_gate, build_validation_policy_preview

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Allowed small-molecule losses by reaction category
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Small-molecule loss tables
# ---------------------------------------------------------------------------
#
# Two-tier design:
#   Tier 1 — _SMALL_MOL_LOSSES: category-specific allow-list (fast path).
#            When the reaction category is known, only these losses are tried.
#   Tier 2 — _UNIVERSAL_LOSSES: category-agnostic fallback.
#            When Tier 1 doesn't fully explain the imbalance (or no category
#            is provided), we greedily try ALL common small-molecule losses.
#            Constraint: universal fallback only explains non-skeleton atoms
#            (H, O, halogen, etc.).  C/N/S differences are never auto-excused.
#
# Each loss entry is a Counter of atoms in the lost molecule.

_LOSS_COUNTERS: Dict[str, Counter] = {
    "H2O":  Counter({"O": 1, "H": 2}),
    "HCl":  Counter({"H": 1, "Cl": 1}),
    "HBr":  Counter({"H": 1, "Br": 1}),
    "HI":   Counter({"H": 1, "I": 1}),
    "HF":   Counter({"H": 1, "F": 1}),
    "H2":   Counter({"H": 2}),
    "N2":   Counter({"N": 2}),
    "CO2":  Counter({"C": 1, "O": 2}),
    "MeOH": Counter({"C": 1, "O": 1, "H": 4}),
    "EtOH": Counter({"C": 2, "O": 1, "H": 6}),
    "AcOH": Counter({"C": 2, "O": 2, "H": 4}),
    "NH3":  Counter({"N": 1, "H": 3}),
    "SO2":  Counter({"S": 1, "O": 2}),
}

# Tier 1: category → allowed loss keys (precise, no false positives)
_SMALL_MOL_LOSSES: Dict[str, List[str]] = {
    "ester_hydrolysis": ["H2O"],
    "amide_formation": ["H2O"],
    "condensation": ["H2O"],
    "aldol": ["H2O"],
    "claisen": ["H2O"],
    "acylation": ["H2O", "HCl"],
    "friedel_crafts": ["H2O", "HCl", "H2"],
    "halogenation": ["HCl", "HBr", "HI"],
    "heck": ["HCl", "HBr", "HI"],
    "elimination": ["H2O", "HCl", "HBr"],
    "diazotization": ["N2"],
    "lactonization": ["H2O"],
    "lactam": ["H2O"],
    "cyclization": ["H2O", "H2"],
    "intramolecular": ["H2O", "HCl", "H2"],
    "ring_formation": ["H2O", "H2"],
    "ring_closing": ["H2O", "H2"],
    "decarboxylation": ["CO2"],
    "curtius": ["CO2", "N2"],
}

# Tier 2: universal fallback — all losses that DON'T involve skeleton atoms.
# These are safe to try without knowing the reaction type.
_UNIVERSAL_LOSS_KEYS: List[str] = [
    "H2O", "HCl", "HBr", "HI", "HF", "H2",
]

_SKELETON_ATOMS: frozenset = frozenset({"C", "N", "S"})

# MCS level mapping (reaction category -> level)
_MCS_REACTION_LEVELS: Dict[str, int] = {
    "ester": 1, "amide": 1, "ether": 1, "hydrolysis": 1,
    "deprotection": 1, "cleavage": 1,
    "alkylation": 2, "acylation": 2, "n_alkylation": 2, "o_alkylation": 2,
    "grignard": 2, "nucleophilic_addition": 2, "electrophilic_addition": 2,
    "halogenation": 2, "oxidation": 2, "reduction": 2,
    "coupling": 2, "cross_coupling": 2, "pd_catalysis": 2,
    "substitution": 2, "elimination": 2,
    "ring_formation": 3, "cyclization": 3, "ring_closing": 3,
    "robinson": 3, "diels_alder": 3, "rcm": 3, "metathesis": 3,
    "aldol": 3, "claisen": 3, "michael": 3,
    "annulation": 3, "cycloaddition": 3, "intramolecular": 3,
    "ring_opening": 3, "ring_retro": 3, "fused_ring": 3,
    "nazarov": 3, "pauson_khand": 3, "fischer_indole": 3,
    "heck_cyclization": 3, "beckmann": 3, "ring_reaction": 3,
    "lactonization": 3, "lactam": 3, "epoxide": 3,
    "baeyer_villiger": 3, "simmons_smith": 3, "cyclopropan": 3,
    "dipolar": 3, "1_3_dipolar": 3, "hantzsch": 3,
    "friedel_crafts_cyclization": 3, "paal_knorr": 3, "knorr_pyrrole": 3,
}

_MCS_LEVEL_CONFIG: Dict[int, Dict[str, Any]] = {
    1: {"coverage_threshold": 0.40, "min_mcs_atoms": 2, "small_mol_exempt": 4},
    2: {"coverage_threshold": 0.30, "min_mcs_atoms": 2, "small_mol_exempt": 6},
    3: {"coverage_threshold": 0.20, "min_mcs_atoms": 1, "small_mol_exempt": 8},
}

# Category -> reaction families for FG compatibility matrix lookup
_CATEGORY_FAMILY_MAP: Dict[str, List[str]] = {
    "oxidation": ["strong_oxidation", "mild_oxidation"],
    "reduction": ["strong_reduction", "mild_reduction"],
    "halogenation": ["halogenation"],
    "cross_coupling": ["pd_catalysis"],
    "coupling": ["pd_catalysis"],
    "suzuki": ["pd_catalysis"],
    "heck": ["pd_catalysis"],
    "grignard": ["strong_base", "organometallic"],
    "alkylation": ["strong_base"],
    "elimination": ["strong_base"],
    "ester_hydrolysis": ["acidic_conditions", "basic_conditions"],
    "amide_formation": ["acidic_conditions"],
    "condensation": ["acidic_conditions"],
    "diazotization": ["acidic_conditions"],
    "friedel_crafts": ["lewis_acid", "electrophilic_aromatic", "strong_acid"],
    "friedel_crafts_alkylation": ["lewis_acid", "electrophilic_aromatic", "strong_acid"],
    "friedel_crafts_acylation": ["lewis_acid", "electrophilic_aromatic", "strong_acid"],
    "nitration": ["strong_acid", "electrophilic_aromatic"],
    "sulfonation": ["strong_acid", "electrophilic_aromatic"],
}



# ---------------------------------------------------------------------------
# _sanitize_product  (three-level fallback, same pattern as bond_break.py)
# ---------------------------------------------------------------------------

def _sanitize_product(mol: Chem.Mol) -> Tuple[bool, Optional[str]]:
    """Apply three-level SanitizeMol fallback strategy.

    Returns ``(success, error_message)``.
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
# check_atom_balance  (public)
# ---------------------------------------------------------------------------

def _get_allowed_losses(reaction_category: Optional[str]) -> List[Counter]:
    """Return Tier-1 allowed loss Counters for the category.

    Returns None if no category-specific entry is found (caller should
    fall back to Tier-2 universal losses).
    """
    if not reaction_category:
        return []
    cat = reaction_category.lower().replace("-", "_").replace(" ", "_")
    keys = _SMALL_MOL_LOSSES.get(cat)
    if keys is None:
        # Try substring match
        for k, v in _SMALL_MOL_LOSSES.items():
            if k in cat or cat in k:
                keys = v
                break
    if keys is None:
        return []
    return [_LOSS_COUNTERS[k] for k in keys if k in _LOSS_COUNTERS]


def _get_universal_losses() -> List[Counter]:
    """Return Tier-2 universal loss Counters (non-skeleton only)."""
    return [_LOSS_COUNTERS[k] for k in _UNIVERSAL_LOSS_KEYS if k in _LOSS_COUNTERS]


def _subtract_losses(
    remaining: Counter,
    loss_counters: List[Counter],
) -> Counter:
    """Greedily subtract as many copies of each loss molecule as possible.

    Modifies and returns *remaining* in-place.  Each loss is tried
    independently; order doesn't matter because they don't share elements
    in practice (H₂O vs HCl vs H₂ etc.).
    """
    for loss_counter in loss_counters:
        if not loss_counter:
            continue
        # How many whole copies of this loss can we subtract?
        max_copies = None
        for elem, cnt in loss_counter.items():
            if cnt <= 0:
                continue
            available = remaining.get(elem, 0)
            copies = available // cnt
            if max_copies is None or copies < max_copies:
                max_copies = copies
        if max_copies and max_copies > 0:
            for elem, cnt in loss_counter.items():
                remaining[elem] -= cnt * max_copies
                if remaining[elem] <= 0:
                    del remaining[elem]
    return remaining


def check_atom_balance(
    precursors: List[str],
    target: str,
    byproduct_smiles: Optional[List[str]] = None,
    reaction_category: Optional[str] = None,
) -> Dict[str, Any]:
    """Atom conservation check.

    Product side = target + byproducts.
    Allows small-molecule losses (H2O/HCl/HBr/N2 etc.) based on
    *reaction_category*.  Skeleton atoms (C/N/S) not conserved triggers
    ``skeleton_imbalance``.  Non-H atom difference > 4 triggers
    ``severe_imbalance``.

    Returns
    -------
    dict
        ``{"balanced", "precursor_atoms", "product_atoms", "excess",
        "deficit", "severe_imbalance", "skeleton_imbalance", "note"}``
    """
    # --- Parse precursors ---
    precursor_atoms: Counter = Counter()
    for smi in precursors:
        mol = parse_mol(smi)
        if mol is None:
            return {
                "ok": False,
                "error": f"invalid precursor SMILES: {smi}",
                "input": smi,
            }
        precursor_atoms += mol_formula_counter(mol)

    # --- Parse target ---
    target_mol = parse_mol(target)
    if target_mol is None:
        return {"ok": False, "error": "invalid target SMILES", "input": target}

    product_atoms: Counter = Counter()
    product_atoms += mol_formula_counter(target_mol)

    # --- Parse and add byproducts ---
    if byproduct_smiles:
        for bp_smi in byproduct_smiles:
            bp_mol = parse_mol(bp_smi.strip())
            if bp_mol is not None:
                product_atoms += mol_formula_counter(bp_mol)
            # silently skip invalid byproduct SMILES

    # --- Compute raw difference ---
    all_elements = set(precursor_atoms.keys()) | set(product_atoms.keys())
    excess: Counter = Counter()   # product has more
    deficit: Counter = Counter()  # precursor has more

    for elem in all_elements:
        diff = product_atoms[elem] - precursor_atoms[elem]
        if diff > 0:
            excess[elem] = diff
        elif diff < 0:
            deficit[elem] = -diff

    # --- Try to explain imbalance via small-molecule losses ---
    #
    # Tier 1: category-specific losses (precise)
    # Tier 2: universal fallback (H2O, HX, H2 — non-skeleton only)
    #
    # Both tiers are applied to both deficit and excess sides.
    remaining_excess = Counter(excess)
    remaining_deficit = Counter(deficit)

    tier1 = _get_allowed_losses(reaction_category)
    universal = _get_universal_losses()

    # Apply Tier 1 first (if available), then Tier 2 for anything left over.
    # Deficit = precursor has more atoms → lost as small molecules
    _subtract_losses(remaining_deficit, tier1)
    _subtract_losses(remaining_deficit, universal)

    # Excess = product has more atoms → gained from small molecules (e.g. hydrolysis)
    _subtract_losses(remaining_excess, tier1)
    _subtract_losses(remaining_excess, universal)

    # --- Determine balance status ---
    # Use adjusted remainders: if all imbalance is explained by allowed
    # small-molecule losses, the reaction is considered balanced.
    balanced = (len(remaining_excess) == 0 and len(remaining_deficit) == 0)

    # Skeleton imbalance check (uses RAW diff).
    #
    # Key insight: in real reactions, precursors often carry "extra" atoms
    # that leave as byproducts (Ph₃P=O in Wittig, Bu₃SnCl in Stille, etc.).
    # This means precursor-side excess of C/N/S is *normal* and expected.
    #
    # What's NOT normal is the product having MORE skeleton atoms than the
    # precursors — that means atoms appeared from nowhere.
    #
    # So we only flag skeleton_imbalance when the PRODUCT side has excess
    # C, N, or S (i.e. the `excess` counter, which = product - precursor > 0).
    skeleton_imbalance = False
    for elem in _SKELETON_ATOMS:
        if excess.get(elem, 0) > 0:
            skeleton_imbalance = True
            break

    # Severe imbalance: total non-H atom difference > 4 (raw).
    # But only count product-side excess (same reasoning: precursor excess
    # is normal for reagent-mediated reactions).  Precursor-side excess of
    # non-H atoms is penalized more softly via the continuous score in
    # validate_forward.
    non_h_product_excess = 0
    for elem in all_elements:
        if elem == "H":
            continue
        non_h_product_excess += excess.get(elem, 0)
    severe_imbalance = non_h_product_excess > 4

    # --- Continuous balance score ---
    # Instead of binary balanced/not, compute a 0-1 score reflecting how
    # much of the imbalance remains unexplained after loss subtraction.
    # This feeds into validate_forward's weighted average.
    total_heavy = sum(v for e, v in precursor_atoms.items() if e != "H")
    total_heavy = max(total_heavy, 1)  # avoid div-by-zero
    unexplained = sum(remaining_deficit.values()) + sum(remaining_excess.values())
    # Ratio of unexplained atoms to total heavy atoms, clamped to [0, 1]
    unexplained_ratio = min(unexplained / total_heavy, 1.0)
    balance_score = round(1.0 - unexplained_ratio, 4)

    # --- Build note ---
    note_parts: List[str] = []
    if balanced:
        if len(excess) == 0 and len(deficit) == 0:
            note_parts.append("atoms fully conserved")
        else:
            note_parts.append("balanced after accounting for small-molecule loss")
    else:
        if remaining_excess:
            note_parts.append(
                "adjusted excess: "
                + ", ".join(f"{e}={c}" for e, c in sorted(remaining_excess.items()))
            )
        if remaining_deficit:
            note_parts.append(
                "adjusted deficit: "
                + ", ".join(f"{e}={c}" for e, c in sorted(remaining_deficit.items()))
            )
        if skeleton_imbalance:
            note_parts.append("product has more skeleton atoms (C/N/S) than precursors")
        if severe_imbalance:
            note_parts.append("severe: product non-H excess > 4")

    return {
        "balanced": balanced,
        "balance_score": balance_score,
        "precursor_atoms": precursor_atoms,
        "product_atoms": product_atoms,
        "excess": excess,
        "deficit": deficit,
        "adjusted_excess": dict(remaining_excess) if remaining_excess else {},
        "adjusted_deficit": dict(remaining_deficit) if remaining_deficit else {},
        "severe_imbalance": severe_imbalance,
        "skeleton_imbalance": skeleton_imbalance,
        "note": "; ".join(note_parts) if note_parts else "",
    }



# ---------------------------------------------------------------------------
# _execute_forward_template
# ---------------------------------------------------------------------------

def _execute_forward_template(
    precursors: List[str],
    target: str,
    template_id: Optional[str] = None,
    reaction_smarts: Optional[str] = None,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """Attempt forward template execution, compare products to target.

    Returns dict with ``attempted``, ``target_in_products``,
    ``tanimoto_to_target``, ``products_generated``, ``best_match``.
    """
    result: Dict[str, Any] = {"attempted": False}

    precursor_mols = [parse_mol(s) for s in precursors]
    target_mol = parse_mol(target)
    if target_mol is None or any(m is None for m in precursor_mols):
        result["error"] = "invalid SMILES in precursors or target"
        return result

    target_can = canonical(target)

    # --- Resolve rxn SMARTS ---
    rxn_smarts: Optional[str] = None

    if template_id:
        try:
            reactions = load_template(
                "reactions.json", knowledge_profile=knowledge_profile
            )
        except Exception:
            reactions = {}
        tmpl = reactions.get(template_id)
        if tmpl and isinstance(tmpl, dict):
            rxn_smarts = tmpl.get("rxn")
            result["template_id"] = template_id
            # If retro template, reverse it for forward execution
            if tmpl.get("type") == "retro" and rxn_smarts and ">>" in rxn_smarts:
                lhs, rhs = rxn_smarts.split(">>", 1)
                rxn_smarts = f"{rhs}>>{lhs}"
                result["smarts_reversed"] = True

    if not rxn_smarts and reaction_smarts:
        rxn_smarts = reaction_smarts

    if not rxn_smarts:
        # No template available — just compute Tanimoto between precursors and target
        best_sim = 0.0
        for pmol in precursor_mols:
            if pmol is not None:
                sim = tanimoto(pmol, target_mol)
                if sim > best_sim:
                    best_sim = sim
        result["tanimoto_to_target"] = round(best_sim, 4)
        result["target_in_products"] = False
        result["products_generated"] = []
        result["best_match"] = ""
        return result

    # --- Execute forward reaction ---
    result["attempted"] = True
    try:
        rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
        n_templates = rxn.GetNumReactantTemplates()
    except Exception as e:
        result["error"] = f"reaction SMARTS parse failed: {e}"
        result["target_in_products"] = False
        result["tanimoto_to_target"] = 0.0
        result["products_generated"] = []
        result["best_match"] = ""
        return result

    # Suppress RDKit warnings during reaction execution
    _lg = RDLogger.logger()
    _lg.setLevel(RDLogger.CRITICAL)
    try:
        raw_products: list = []
        valid_mols = [m for m in precursor_mols if m is not None]

        if n_templates == 1:
            for pmol in valid_mols:
                try:
                    prods = rxn.RunReactants((pmol,))
                    raw_products.extend(prods)
                except Exception:
                    pass
        elif n_templates == 2 and len(valid_mols) >= 2:
            for i in range(len(valid_mols)):
                for j in range(len(valid_mols)):
                    if i == j:
                        continue
                    try:
                        prods = rxn.RunReactants((valid_mols[i], valid_mols[j]))
                        raw_products.extend(prods)
                    except Exception:
                        pass
        elif n_templates == 2 and len(valid_mols) == 1:
            try:
                prods = rxn.RunReactants((valid_mols[0], valid_mols[0]))
                raw_products.extend(prods)
            except Exception:
                pass
    finally:
        _lg.setLevel(RDLogger.ERROR)

    # Sanitize and collect unique products
    seen: Set[str] = set()
    valid_products: List[str] = []
    for prod_tuple in raw_products:
        try:
            all_ok = all(_sanitize_product(p)[0] for p in prod_tuple)
            if not all_ok:
                continue
            smi = ".".join(Chem.MolToSmiles(p) for p in prod_tuple)
            if smi and smi not in seen:
                seen.add(smi)
                valid_products.append(smi)
        except Exception:
            pass

    # Compare products to target
    target_in_products = False
    best_tanimoto = 0.0
    best_match = ""

    for prod_smi in valid_products[:20]:
        for frag_smi in prod_smi.split("."):
            frag_can = canonical(frag_smi)
            if frag_can == target_can:
                target_in_products = True
                best_tanimoto = 1.0
                best_match = frag_can or frag_smi
                break
            frag_mol = parse_mol(frag_smi)
            if frag_mol is not None and target_mol is not None:
                sim = tanimoto(frag_mol, target_mol)
                if sim > best_tanimoto:
                    best_tanimoto = sim
                    best_match = frag_smi
        if target_in_products:
            break

    result["target_in_products"] = target_in_products
    result["tanimoto_to_target"] = round(best_tanimoto, 4)
    result["products_generated"] = valid_products[:10]
    result["best_match"] = best_match
    return result



# ---------------------------------------------------------------------------
# _check_scaffold_alignment  (MCS-based)
# ---------------------------------------------------------------------------

def _get_mcs_level(reaction_category: Optional[str]) -> int:
    """Map reaction category to MCS check level (1=strict, 2=medium, 3=loose)."""
    if not reaction_category:
        return 1
    cat = reaction_category.lower().replace("-", "_").replace(" ", "_")
    if cat in _MCS_REACTION_LEVELS:
        return _MCS_REACTION_LEVELS[cat]
    for keyword, level in _MCS_REACTION_LEVELS.items():
        if keyword in cat or cat in keyword:
            return level
    return 1


def _check_scaffold_alignment(
    precursors: List[str],
    target: str,
    reaction_category: Optional[str] = None,
) -> Dict[str, Any]:
    """MCS-based scaffold alignment check.

    Returns dict with ``aligned``, ``coverage_ratio``, ``mcs_level``,
    ``warning``.
    """
    try:
        from rdkit.Chem import rdFMCS
    except ImportError:
        return {
            "aligned": True,
            "coverage_ratio": 0.0,
            "mcs_level": 1,
            "warning": "rdFMCS not available",
        }

    target_mol = parse_mol(target)
    if target_mol is None:
        return {
            "aligned": False,
            "coverage_ratio": 0.0,
            "mcs_level": 1,
            "warning": "invalid target SMILES",
        }

    target_heavy = target_mol.GetNumHeavyAtoms()
    if target_heavy == 0:
        return {"aligned": True, "coverage_ratio": 1.0, "mcs_level": 1, "warning": None}

    mcs_level = _get_mcs_level(reaction_category)
    config = _MCS_LEVEL_CONFIG[mcs_level]

    total_mcs_atoms = 0
    warning = None

    for smi in precursors:
        pmol = parse_mol(smi)
        if pmol is None:
            continue
        p_heavy = pmol.GetNumHeavyAtoms()

        try:
            if mcs_level >= 3:
                mcs_result = rdFMCS.FindMCS(
                    [pmol, target_mol],
                    timeout=5,
                    matchValences=False,
                    ringMatchesRingOnly=False,
                    completeRingsOnly=False,
                    bondCompare=rdFMCS.BondCompare.CompareAny,
                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                )
            else:
                mcs_result = rdFMCS.FindMCS(
                    [pmol, target_mol],
                    timeout=5,
                    matchValences=False,
                    ringMatchesRingOnly=True,
                    completeRingsOnly=False,
                    bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                )
            mcs_atoms = mcs_result.numAtoms if mcs_result else 0
        except Exception:
            mcs_atoms = 0

        # Check minimum MCS atoms for non-small molecules
        if p_heavy >= config["small_mol_exempt"] and mcs_atoms < config["min_mcs_atoms"]:
            warning = (
                f"precursor {smi} (heavy={p_heavy}) has MCS of only "
                f"{mcs_atoms} atoms with target (Level {mcs_level} "
                f"requires >= {config['min_mcs_atoms']})"
            )

        total_mcs_atoms += mcs_atoms

    coverage = min(total_mcs_atoms / target_heavy, 1.0)
    aligned = coverage >= config["coverage_threshold"]

    if not aligned and warning is None:
        warning = (
            f"MCS coverage {coverage:.1%} below threshold "
            f"{config['coverage_threshold']:.0%} (Level {mcs_level})"
        )

    return {
        "aligned": aligned,
        "coverage_ratio": round(coverage, 3),
        "mcs_level": mcs_level,
        "warning": warning,
    }



# ---------------------------------------------------------------------------
# _check_bond_topology
# ---------------------------------------------------------------------------

# Polar addition / condensation categories that cannot form bonds on
# aromatic atoms.
_POLAR_ADDITION_CATS: frozenset = frozenset({
    "aldol", "claisen", "michael", "mannich", "wittig",
    "grignard", "nucleophilic_addition", "condensation",
    "amide_formation", "ester", "esterification",
    "robinson", "knoevenagel", "henry", "reformatsky",
})

# Categories that CAN form aromatic bonds (cross-coupling, EAS, etc.)
_AROMATIC_BOND_CATS: frozenset = frozenset({
    "cross_coupling", "coupling", "suzuki", "heck", "sonogashira",
    "buchwald", "negishi", "stille", "kumada",
    "friedel_crafts", "electrophilic_aromatic_substitution",
    "c_h_activation", "fischer_indole", "hantzsch",
    "paal_knorr", "knorr_pyrrole",
})

_PAAL_KNORR_CATS: frozenset = frozenset({
    "paal_knorr",
    "knorr_pyrrole",
    "paal_knorr_pyrrole",
})


def _normalize_reaction_category(reaction_category: Optional[str]) -> str:
    if not reaction_category:
        return ""
    return reaction_category.lower().replace("-", "_").replace(" ", "_")


def _category_matches(cat: str, key: str) -> bool:
    """Match reaction-category keys on token boundaries, not raw substrings.

    Raw substring matching makes short keys such as ``ester`` hit unrelated
    words like ``thioester``.  That is too aggressive for hard-block checks.
    """
    if not cat or not key:
        return False
    if cat == key:
        return True

    cat_tokens = [token for token in cat.split("_") if token]
    key_tokens = [token for token in key.split("_") if token]
    if not cat_tokens or not key_tokens:
        return False
    if len(key_tokens) == 1:
        return key_tokens[0] in cat_tokens

    window = len(key_tokens)
    return any(cat_tokens[i : i + window] == key_tokens for i in range(len(cat_tokens) - window + 1))


def _category_matches_any(cat: str, keys: Iterable[str]) -> bool:
    return any(_category_matches(cat, key) for key in keys)


def _count_aromatic_rings(mol: Chem.Mol) -> int:
    count = 0
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            count += 1
    return count


def _has_fused_aromatic_system(mol: Chem.Mol) -> bool:
    aromatic_rings: List[Set[int]] = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        ring_set = set(ring)
        if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            aromatic_rings.append(ring_set)

    for i in range(len(aromatic_rings)):
        for j in range(i + 1, len(aromatic_rings)):
            if len(aromatic_rings[i] & aromatic_rings[j]) >= 2:
                return True
    return False


def _custom_or_declared_topology(action_context: Optional[Dict[str, Any]]) -> bool:
    action_context = action_context or {}
    risk_tags = action_context.get("risk_tags") or []
    return bool(
        action_context.get("source") == "custom_precursors"
        or action_context.get("changed_bonds")
        or action_context.get("preserved_anchors")
        or action_context.get("family_evidence")
        or "route_sketch_derived" in risk_tags
    )


def _check_reaction_specific_plausibility(
    precursors: List[str],
    target: str,
    reaction_category: Optional[str] = None,
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Reaction-family specific plausibility checks.

    Current rule set focuses on preventing unrealistic Paal-Knorr overreach:
    Paal-Knorr / Knorr pyrrole reactions should not create a fused multi-aromatic
    scaffold from non-fused / non-aromatic precursors in a single step.
    """
    cat = _normalize_reaction_category(reaction_category)
    if not cat:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "no reaction category, skipped",
        }

    is_paal_knorr = _category_matches_any(cat, _PAAL_KNORR_CATS)
    if not is_paal_knorr:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "no reaction-specific restrictions triggered",
        }

    target_mol = parse_mol(target)
    if target_mol is None:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "invalid target, skipped",
        }

    violations: List[Dict[str, Any]] = []
    custom_or_declared_topology = _custom_or_declared_topology(action_context)
    violation_severity = "requires_evidence" if custom_or_declared_topology else "hard_fail"
    target_arom_rings = _count_aromatic_rings(target_mol)
    target_fused_arom = _has_fused_aromatic_system(target_mol)

    precursor_arom_ring_counts: List[int] = []
    precursor_has_fused_arom = False
    for smi in precursors:
        pmol = parse_mol(smi)
        if pmol is None:
            continue
        precursor_arom_ring_counts.append(_count_aromatic_rings(pmol))
        if _has_fused_aromatic_system(pmol):
            precursor_has_fused_arom = True

    max_prec_arom = max(precursor_arom_ring_counts) if precursor_arom_ring_counts else 0
    sum_prec_arom = sum(precursor_arom_ring_counts)

    # Hard fail: creating fused aromatic scaffolds from non-fused/non-aromatic inputs
    if target_fused_arom and (not precursor_has_fused_arom) and max_prec_arom <= 1 and sum_prec_arom <= 1:
        violations.append({
            "reason": (
                f"reaction '{cat}' implies Paal-Knorr pyrrole formation, but product has fused aromatic scaffold "
                "while precursors contain no fused aromatic system"
            ),
            "severity": "hard_fail",
            "type": "paal_knorr_fused_aromatic_overreach",
        })

    # Hard fail: generating 2+ aromatic rings in one Paal-Knorr step from fully non-aromatic precursors
    if target_arom_rings >= 2 and sum_prec_arom == 0:
        violations.append({
            "reason": (
                f"reaction '{cat}' would create {target_arom_rings} aromatic ring(s) from non-aromatic precursors "
                "in a single Paal-Knorr step"
            ),
            "severity": "hard_fail",
            "type": "paal_knorr_aromatic_ring_overcreation",
        })

    has_hard_fail = any(v.get("severity") == "hard_fail" for v in violations)
    return {
        "pass": not has_hard_fail,
        "has_hard_fail": has_hard_fail,
        "violations": violations,
        "summary": "ok" if not violations else f"{len(violations)} reaction-specific violation(s)",
    }


def _check_bond_topology(
    precursors: List[str],
    target: str,
    reaction_category: Optional[str] = None,
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Bond change topology check.

    Compares ring structures between precursors and target to detect
    infeasible aromatic ring formation for polar addition reactions.

    Returns dict with ``pass``, ``has_hard_fail``, ``violations``,
    ``summary``.
    """
    action_context = action_context or {}
    target_mol = parse_mol(target)
    if target_mol is None:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "invalid target, skipped",
        }

    if not reaction_category:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "no reaction category, skipped",
        }

    cat = reaction_category.lower().replace("-", "_").replace(" ", "_")
    custom_or_declared_topology = _custom_or_declared_topology(action_context)
    violation_severity = (
        "requires_evidence"
        if custom_or_declared_topology
        else "hard_fail"
    )

    # Only check polar addition reactions that cannot form aromatic bonds
    is_polar = _category_matches_any(cat, _POLAR_ADDITION_CATS)
    can_aromatic = _category_matches_any(cat, _AROMATIC_BOND_CATS)

    if not is_polar or can_aromatic:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "bond topology consistent with reaction type",
        }

    violations: List[Dict[str, Any]] = []

    # Check: new aromatic rings in product not present in precursors
    precursor_aromatic_count = 0
    for smi in precursors:
        pmol = parse_mol(smi)
        if pmol is None:
            continue
        ri = pmol.GetRingInfo()
        for ring in ri.AtomRings():
            if all(pmol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
                precursor_aromatic_count += 1

    target_ri = target_mol.GetRingInfo()
    target_aromatic_count = 0
    for ring in target_ri.AtomRings():
        if all(target_mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            target_aromatic_count += 1

    new_aromatic = target_aromatic_count - precursor_aromatic_count
    if new_aromatic > 0:
        if custom_or_declared_topology:
            reason = (
                f"product has {new_aromatic} more aromatic ring(s) than precursors; "
                f"reaction label '{cat}' does not by itself establish this aromatic "
                "ring construction and needs atom-source/aromatization proof"
            )
        else:
            reason = (
                f"polar addition/condensation reaction '{cat}' cannot form "
                f"aromatic rings, but product has {new_aromatic} more "
                f"aromatic ring(s) than precursors"
            )
        violations.append({
            "reason": reason,
            "severity": violation_severity,
            "type": "aromatic_ring_formation",
            "code": "aromatic_ring_formation",
        })

    # Check: aromatic-nonaromatic ring fusion in product not in precursors
    target_arom_atoms: Set[int] = set()
    target_non_arom_ring_atoms: Set[int] = set()
    for ring in target_ri.AtomRings():
        atoms_set = set(ring)
        if all(target_mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            target_arom_atoms |= atoms_set
        else:
            target_non_arom_ring_atoms |= atoms_set

    fusion_atoms = target_arom_atoms & target_non_arom_ring_atoms
    if len(fusion_atoms) >= 2:
        # Check if precursors already have this fusion
        precursor_has_fusion = False
        for smi in precursors:
            pmol = parse_mol(smi)
            if pmol is None:
                continue
            pri = pmol.GetRingInfo()
            p_arom: Set[int] = set()
            p_non_arom: Set[int] = set()
            for ring in pri.AtomRings():
                ring_set = set(ring)
                if all(pmol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
                    p_arom |= ring_set
                else:
                    p_non_arom |= ring_set
            if p_arom & p_non_arom:
                precursor_has_fusion = True
                break

        if not precursor_has_fusion:
            if custom_or_declared_topology:
                reason = (
                    f"product has aromatic-nonaromatic ring fusion "
                    f"(fusion atoms: {sorted(fusion_atoms)}) not present "
                    f"in precursors; reaction label '{cat}' does not by itself "
                    "establish the fused scaffold closure and needs atom-source proof "
                    "for the new ring bond and preserved junction atoms"
                )
            else:
                reason = (
                    f"product has aromatic-nonaromatic ring fusion "
                    f"(fusion atoms: {sorted(fusion_atoms)}) not present "
                    f"in precursors; polar reaction '{cat}' cannot form "
                    f"bonds on aromatic atoms to build fused scaffold"
                )
            violations.append({
                "reason": reason,
                "severity": violation_severity,
                "type": "fusion_bond_mapping_requires_evidence",
                "code": "fusion_bond_mapping_requires_evidence",
            })

    has_hard_fail = any(v["severity"] == "hard_fail" for v in violations)

    if not violations:
        summary = "bond topology consistent with reaction type"
    else:
        summary = f"{len(violations)} bond topology violation(s)"

    return {
        "pass": not has_hard_fail,
        "has_hard_fail": has_hard_fail,
        "violations": violations,
        "summary": summary,
        "policy_note": (
            "custom topology action requires LLM atom-source proof"
            if violations and custom_or_declared_topology
            else ""
        ),
    }



# ---------------------------------------------------------------------------
# _check_fg_compatibility
# ---------------------------------------------------------------------------

def _check_fg_compatibility(
    precursors: List[str],
    target: str,
    reaction_category: Optional[str] = None,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """Functional group compatibility check using fg_compatibility_matrix.json.

    Returns dict with ``compatible``, ``warnings``.
    """
    if not reaction_category:
        return {"compatible": True, "warnings": []}

    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    warnings: List[Dict[str, Any]] = []
    knowledge_refs: List[Dict[str, Any]] = []
    compatible = True
    for smiles in list(precursors) + [target]:
        result = check_reaction_compatibility(
            smiles,
            reaction_category,
            knowledge_profile=knowledge_profile,
        )
        if "error" in result:
            continue
        compatible = compatible and result.get("compatible", True)
        warnings.extend(result.get("warnings", []))
        knowledge_refs.extend(result.get("knowledge_refs", []) or [])

    # Deduplicate warnings by (fg, level)
    seen: Set[str] = set()
    deduped: List[Dict[str, Any]] = []
    for w in warnings:
        key = f"{w['fg']}|{w['level']}|{w['reaction_family']}"
        if key not in seen:
            seen.add(key)
            deduped.append(w)

    deduped_refs: List[Dict[str, Any]] = []
    seen_refs: Set[Tuple[str, ...]] = set()
    for ref in knowledge_refs:
        if not isinstance(ref, dict):
            continue
        key = tuple(f"{name}={ref[name]}" for name in sorted(ref))
        if key not in seen_refs:
            seen_refs.add(key)
            deduped_refs.append(ref)

    return {
        "compatible": compatible,
        "warnings": deduped,
        "knowledge_profile_hash": knowledge_profile.digest,
        "knowledge_refs": deduped_refs,
    }


def _audit_precursor_state(
    precursors: List[str],
    action_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Report open-shell precursor atoms without inferring chemistry from names."""
    rows: List[Dict[str, Any]] = []
    elemental_metals: List[Dict[str, Any]] = []
    total_radical_electrons = 0
    context = action_context or {}
    reagent_smiles = {
        canonical(smiles) or smiles
        for smiles in (context.get("reagents") or [])
    }
    for smiles in precursors:
        mol = parse_mol(smiles)
        if mol is None:
            continue
        if mol.GetNumAtoms() == 1:
            atom = mol.GetAtomWithIdx(0)
            if atom.GetSymbol() in {"Li", "Mg", "Zn", "Cu"}:
                normalized = canonical(smiles) or smiles
                elemental_metals.append({
                    "smiles": normalized,
                    "element": atom.GetSymbol(),
                    "role": "reagent" if normalized in reagent_smiles else "precursor",
                    "representation_note": (
                        "single-atom elemental metal encoding; radical electrons are "
                        "an RDKit representation fact, not an unsupported molecular radical"
                    ),
                })
                continue
        atoms = []
        for atom in mol.GetAtoms():
            radical_electrons = int(atom.GetNumRadicalElectrons())
            if radical_electrons <= 0:
                continue
            total_radical_electrons += radical_electrons
            atoms.append({
                "atom_idx": atom.GetIdx(),
                "element": atom.GetSymbol(),
                "radical_electrons": radical_electrons,
                "formal_charge": atom.GetFormalCharge(),
            })
        if atoms:
            rows.append({
                "smiles": canonical(smiles) or smiles,
                "radical_electrons": sum(
                    atom["radical_electrons"] for atom in atoms
                ),
                "atoms": atoms,
            })

    structured_tokens = {
        str(item or "").strip().lower()
        for key in ("risk_tags", "intended_deltas")
        for item in (context.get(key) or [])
    }
    family_evidence = context.get("family_evidence") or {}
    declared = bool(
        structured_tokens
        & {
            "open_shell_precursor",
            "radical_precursor",
            "radical_intermediate",
            "radical_capture",
        }
        or (
            isinstance(family_evidence, dict)
            and any(
                key in family_evidence
                for key in (
                    "open_shell_precursor_role",
                    "radical_precursor_role",
                    "radical_atom_source",
                )
            )
        )
    )
    if not rows:
        return {
            "status": "elemental_metal_reagent" if elemental_metals else "closed_shell",
            "declared": declared,
            "total_radical_electrons": 0,
            "open_shell_precursors": [],
            "elemental_metal_reagents": elemental_metals,
        }
    return {
        "status": "open_shell_detected",
        "declared": declared,
        "total_radical_electrons": total_radical_electrons,
        "open_shell_precursors": rows,
        "elemental_metal_reagents": elemental_metals,
    }


# ---------------------------------------------------------------------------
# validate_forward  (public main entry)
# ---------------------------------------------------------------------------

def validate_forward(
    precursors: List[str],
    target: str,
    template_id: Optional[str] = None,
    reaction_smarts: Optional[str] = None,
    reaction_category: Optional[str] = None,
    byproduct_smiles: Optional[List[str]] = None,
    action_context: Optional[Dict[str, Any]] = None,
    reagents: Optional[List[str]] = None,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """Forward synthesis feasibility validation (main entry).

    Runs 13 sequential checks:
      1. Atom balance  (check_atom_balance)
      2. Template forward execution  (_execute_forward_template)
      3. MCS scaffold alignment  (_check_scaffold_alignment)
      4. Graph delta  (audit_graph_delta)
      5. Ring topology  (audit_ring_topology)
      6. Bond change topology  (_check_bond_topology)
      7. Atom-source audit  (audit_atom_mapping)
      8. Reaction-specific plausibility (_check_reaction_specific_plausibility)
      9. Functional-group delta  (audit_fg_delta)
      10. Reaction-family interpretation  (validate_reaction_family)
      11. Functional group compatibility  (_check_fg_compatibility)
      12. Precursor electronic state  (_audit_precursor_state)
      13. EAS site-selectivity evidence  (audit_eas_site_selectivity)

    Hard Gate — any one of these triggers ``pass=False``:
      - severe_imbalance from atom balance
      - skeleton_imbalance from atom balance
      - scaffold not aligned (aligned=False)
      - bond topology has_hard_fail
      - reaction-specific plausibility has_hard_fail
      - reaction-family topology has_hard_fail
      - ring-topology hard-fail finding
      - forbidden FG (compatible=False)

    Returns
    -------
    dict
        ``{"ok", "target", "precursors", "checks": {...},
        "assessment": {"feasibility_score", "pass", "hard_fail_reasons",
        "score_breakdown"}}``
    """
    # --- Validate inputs ---
    target_mol = parse_mol(target)
    if target_mol is None:
        return {"ok": False, "error": "invalid target SMILES", "input": target}

    reagents = list(reagents or [])
    for smi in list(precursors) + reagents:
        if parse_mol(smi) is None:
            return {"ok": False, "error": f"invalid reaction-input SMILES: {smi}", "input": smi}

    reaction_inputs = list(precursors) + reagents

    # --- Step 1: Atom balance ---
    atom_bal = check_atom_balance(
        reaction_inputs, target,
        byproduct_smiles=byproduct_smiles,
        reaction_category=reaction_category,
    )
    # check_atom_balance may return {"ok": False, ...} on invalid SMILES
    if atom_bal.get("ok") is False:
        return atom_bal

    # --- Step 2: Template forward execution ---
    template_exec = _execute_forward_template(
        precursors, target,
        template_id=template_id,
        reaction_smarts=reaction_smarts,
        knowledge_profile=knowledge_profile,
    )

    # --- Step 3: MCS scaffold alignment ---
    scaffold_align = _check_scaffold_alignment(
        precursors, target,
        reaction_category=reaction_category,
    )

    # --- Step 4: Generic graph delta ---
    graph_delta = audit_graph_delta(
        target,
        precursors,
        reaction_category=reaction_category,
        action_context=action_context,
    )

    # --- Step 5: Ring topology delta ---
    ring_topology = audit_ring_topology(
        target,
        precursors,
        reaction_category=reaction_category,
        action_context=action_context,
    )

    # --- Step 6: Bond change topology ---
    bond_topo = _check_bond_topology(
        precursors, target,
        reaction_category=reaction_category,
        action_context=action_context,
    )

    # --- Step 7: MCS-derived atom-source audit ---
    atom_mapping = audit_atom_mapping(
        target,
        precursors,
        reaction_category=reaction_category,
        action_context=action_context,
    )

    # --- Step 8: Reaction-specific plausibility ---
    reaction_specific = _check_reaction_specific_plausibility(
        precursors, target,
        reaction_category=reaction_category,
        action_context=action_context,
    )

    # --- Step 9: Functional-group delta ---
    fg_delta = audit_fg_delta(
        target,
        precursors,
        reaction_category=reaction_category,
        action_context=action_context,
    )

    # --- Step 10: Reaction-family topology plausibility ---
    ring_family = validate_reaction_family(
        precursors,
        target,
        reaction_category=reaction_category,
        ring_topology=ring_topology,
        graph_delta=graph_delta,
        fg_delta=fg_delta,
        template_execution=template_exec,
        action_context=action_context,
        knowledge_profile=knowledge_profile,
    )

    # --- Step 11: Functional group compatibility ---
    fg_compat = _check_fg_compatibility(
        precursors, target,
        reaction_category=reaction_category,
        knowledge_profile=knowledge_profile,
    )

    # --- Step 12: Precursor electronic-state audit ---
    precursor_state = _audit_precursor_state(reaction_inputs, action_context)

    # --- Step 13: EAS site-selectivity evidence ---
    site_audit = audit_site_retention(target, precursors)
    eas_site_selectivity = audit_eas_site_selectivity(
        target,
        precursors,
        reaction_type=reaction_category or "",
        action_context=action_context,
        site_audit=site_audit,
        knowledge_profile=knowledge_profile,
    )

    # --- Hard Gate evaluation ---
    hard_fail_reasons: List[str] = []

    if atom_bal.get("severe_imbalance"):
        hard_fail_reasons.append("severe_imbalance")
    if atom_bal.get("skeleton_imbalance"):
        hard_fail_reasons.append("skeleton_imbalance")
    if not scaffold_align.get("aligned", True):
        hard_fail_reasons.append("scaffold_not_aligned")
    if bond_topo.get("has_hard_fail"):
        hard_fail_reasons.append("bond_topology_violation")
    if reaction_specific.get("has_hard_fail"):
        hard_fail_reasons.append("reaction_specific_violation")
    if ring_family.get("has_hard_fail"):
        hard_fail_reasons.append("ring_family_violation")
    if any(
        item.get("severity") == "hard_fail"
        for item in ring_topology.get("violations", []) or []
    ):
        hard_fail_reasons.append("ring_topology_violation")
    if not fg_compat.get("compatible", True):
        hard_fail_reasons.append("forbidden_fg")
    overall_pass = len(hard_fail_reasons) == 0

    # --- Score breakdown ---
    # atom_balance_score: use continuous score from check_atom_balance
    if atom_bal.get("severe_imbalance") or atom_bal.get("skeleton_imbalance"):
        atom_balance_score = 0.0
    else:
        atom_balance_score = float(atom_bal.get("balance_score", 0.5))

    # template_match_score
    template_match_score = float(template_exec.get("tanimoto_to_target", 0.0))

    # scaffold_alignment_score
    scaffold_alignment_score = float(scaffold_align.get("coverage_ratio", 0.0))

    # bond_topology_score
    bond_topology_score = 0.0 if bond_topo.get("has_hard_fail") else 1.0

    # fg_compatibility_score
    if not fg_compat.get("compatible", True):
        fg_compatibility_score = 0.0
    elif fg_compat.get("warnings"):
        # Has risky warnings but no forbidden
        fg_compatibility_score = 0.5
    else:
        fg_compatibility_score = 1.0

    # Weighted average
    feasibility_score = (
        atom_balance_score * 0.25
        + template_match_score * 0.25
        + scaffold_alignment_score * 0.20
        + bond_topology_score * 0.15
        + fg_compatibility_score * 0.15
    )
    feasibility_score = round(max(0.0, min(1.0, feasibility_score)), 4)

    result = {
        "ok": True,
        "target": canonical(target) or target,
        "precursors": [canonical(s) or s for s in precursors],
        "reagents": [canonical(s) or s for s in reagents],
        "checks": {
            "atom_balance": atom_bal,
            "template_execution": template_exec,
            "scaffold_alignment": scaffold_align,
            "graph_delta": graph_delta,
            "ring_topology": ring_topology,
            "bond_topology": bond_topo,
            "atom_mapping": atom_mapping,
            "reaction_specific": reaction_specific,
            "fg_delta": fg_delta,
            "ring_family": ring_family,
            "fg_compatibility": fg_compat,
            "precursor_state": precursor_state,
            "eas_site_selectivity": eas_site_selectivity,
            "action_context": dict(action_context or {}),
        },
        "assessment": {
            "feasibility_score": feasibility_score,
            "pass": overall_pass,
            "hard_fail_reasons": hard_fail_reasons if hard_fail_reasons else None,
            "score_breakdown": {
                "atom_balance": atom_balance_score,
                "template_match": template_match_score,
                "scaffold_alignment": scaffold_alignment_score,
                "bond_topology": bond_topology_score,
                "fg_compatibility": fg_compatibility_score,
            },
        },
    }
    result["assessment"]["gate"] = build_validation_gate(
        result,
        ring_topology_audit=ring_topology,
        action_context=action_context,
        eas_site_audit=eas_site_selectivity,
    )
    result["assessment"]["policy_preview"] = build_validation_policy_preview(result)
    result["assessment"]["evidence_packet"] = build_validation_evidence_packet(result)
    return result
