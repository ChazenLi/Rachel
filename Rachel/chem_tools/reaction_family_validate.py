"""Reaction-family proof obligations for topology-changing steps.

This module translates a proposed reaction family into validation obligations.
It does not prove synthetic feasibility by itself. Registered family names are
optional hints; declared atom-source, tether, and anchor evidence must remain
visible even when the name is unregistered.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set

from .validation_findings import make_finding
from .topology_intent import action_declares_topology_change


def _family_category(reaction_category: Optional[str]) -> str:
    return (reaction_category or "").lower().replace("-", "_").replace(" ", "_")


def _has_any(text: str, needles: Iterable[str]) -> bool:
    return any(needle in text for needle in needles)


def _has_action_token(text: str, tokens: Iterable[str]) -> bool:
    words = {word for word in text.split("_") if word}
    return bool(words & set(tokens))


def _has_complex_ring_delta(ring_topology: Dict[str, Any]) -> bool:
    delta_codes = set(ring_topology.get("delta_codes", []) or [])
    return bool(
        delta_codes
        & {
            "new_fused_ring_system",
            "new_spiro_center",
            "new_bridged_system",
            "new_macrocycle",
            "new_medium_ring",
            "ring_system_merge",
        }
    )


def _delta_codes(payload: Optional[Dict[str, Any]]) -> Set[str]:
    if not isinstance(payload, dict):
        return set()
    return {str(code) for code in (payload.get("delta_codes", []) or []) if str(code)}


def _compact_family_evidence(raw: Any) -> Dict[str, Any]:
    if not isinstance(raw, dict):
        return {}
    return {
        str(key): value
        for key, value in raw.items()
        if value not in (None, "", [], {})
    }


_EVIDENCE_KEY_ALIASES = {
    "dipole": ("dipole", "dipole_atom_source", "one_three_dipole"),
    "dipolarophile": (
        "dipolarophile",
        "dipolarophile_atom_source",
        "alkene_partner",
    ),
    "regiochemistry_rationale": (
        "regiochemistry_rationale",
        "regioselectivity_rationale",
        "regioselectivity",
    ),
    "stereochemistry_rationale": (
        "stereochemistry_rationale",
        "stereoselectivity_rationale",
        "stereochemistry",
    ),
}


def _family_evidence_summary(
    raw: Any,
    required_keys: Iterable[str],
) -> Dict[str, Any]:
    family_evidence = _compact_family_evidence(raw)
    required = list(dict.fromkeys(str(key) for key in required_keys if str(key)))
    if not required:
        provided_keys = list(family_evidence)
        return {
            "state": "provided_unprofiled" if provided_keys else "not_required",
            "provided_keys": provided_keys,
            "missing_keys": [],
            "provided": dict(family_evidence),
            "required_keys": [],
            "missing_evidence_codes": [],
        }

    provided_keys: List[str] = []
    missing_keys: List[str] = []
    provided: Dict[str, Any] = {}
    for key in required:
        aliases = _EVIDENCE_KEY_ALIASES.get(key, (key,))
        matched = next((alias for alias in aliases if alias in family_evidence), "")
        if matched:
            provided_keys.append(matched)
            provided[matched] = family_evidence[matched]
        else:
            missing_keys.append(key)
    if missing_keys and provided_keys:
        state = "partial"
    elif missing_keys:
        state = "missing"
    else:
        state = "provided"
    return {
        "state": state,
        "provided_keys": provided_keys,
        "missing_keys": missing_keys,
        "provided": provided,
        "required_keys": required,
        "missing_evidence_codes": ["family_evidence_missing"] if missing_keys else [],
    }


_PRESERVE_CODES = [
    "skeleton_imbalance",
    "severe_imbalance",
    "olefination_lacks_intramolecular_tether_proof",
    "olefination_ring_closure_template_miss",
    "rcm_requires_single_diene_precursor",
    "new_fused_medium_ring_requires_evidence",
    "new_fused_ring_requires_evidence",
    "new_spiro_center_requires_evidence",
    "new_bridged_system_requires_evidence",
    "new_macrocycle_requires_evidence",
    "apparent_scaffold_jump",
]


def _knowledge_family_profile(
    cat: str,
    action_context: Optional[Dict[str, Any]] = None,
    knowledge_profile=None,
) -> Optional[Dict[str, Any]]:
    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    category = _family_category(cat)
    action_context = action_context or {}
    evidence_keys = set(
        _compact_family_evidence(action_context.get("family_evidence"))
    )
    declared_deltas = {
        str(item or "").strip().lower()
        for item in action_context.get("intended_deltas", []) or []
        if str(item or "").strip()
    }
    high_risk = {
        "new_fused_ring_system",
        "new_spiro_center",
        "new_bridged_system",
        "new_medium_ring",
        "new_macrocycle",
    }
    declared_high_risk = bool(
        action_declares_topology_change(action_context)
        and (
            declared_deltas & high_risk
            or evidence_keys
            & {
                "junction_source",
                "bridgehead_source",
                "high_risk_topology_proof",
            }
        )
    )

    resource = knowledge_profile.get("chem.reaction_families")
    entries = [
        entry
        for entry in resource.get("entries", []) or []
        if isinstance(entry, dict) and entry.get("id")
    ]
    entries.sort(
        key=lambda entry: (
            -int(entry.get("priority", 0) or 0),
            str(entry.get("id", "")),
        )
    )
    category_tokens = set(re.findall(r"[a-z0-9]+", category))
    for entry in entries:
        match = entry.get("match") or {}
        any_terms = [str(item).lower() for item in match.get("any_terms", []) or []]
        any_evidence = set(match.get("any_evidence_keys", []) or [])
        alternatives = []
        if any_terms:
            alternatives.append(any(term in category for term in any_terms))
        if any_evidence:
            alternatives.append(bool(any_evidence & evidence_keys))
        if alternatives and not any(alternatives):
            continue
        all_terms = [str(item).lower() for item in match.get("all_terms", []) or []]
        if any(term not in category for term in all_terms):
            continue
        groups = match.get("all_term_groups", []) or []
        if any(
            not any(str(term).lower() in category for term in group)
            for group in groups
        ):
            continue
        action_tokens = set(match.get("any_action_tokens", []) or [])
        if action_tokens and not action_tokens & category_tokens:
            continue
        if match.get("declared_high_risk_topology") and not declared_high_risk:
            continue
        profile = {
            "family_key": str(entry.get("family_key", "")),
            "family_class": str(entry.get("family_class", "")),
            "allowed_deltas": set(entry.get("allowed_deltas", []) or []),
            "forbidden_deltas": set(entry.get("forbidden_deltas", []) or []),
            "required_evidence": list(entry.get("required_evidence", []) or []),
            "evidence_keys": list(entry.get("evidence_keys", []) or []),
            "risk_level": str(entry.get("risk_level", "medium")),
            "knowledge_ref": knowledge_profile.source(
                "chem.reaction_families", str(entry["id"])
            ),
        }
        return profile
    return None


def _family_profile(
    cat: str,
    action_context: Optional[Dict[str, Any]] = None,
    knowledge_profile=None,
) -> Optional[Dict[str, Any]]:
    """Compatibility adapter backed only by the resolved family resource."""
    return _knowledge_family_profile(
        cat,
        action_context=action_context,
        knowledge_profile=knowledge_profile,
    )


def _build_family_interpretation(
    *,
    cat: str,
    precursors: Sequence[str],
    ring_topology: Dict[str, Any],
    graph_delta: Optional[Dict[str, Any]],
    fg_delta: Optional[Dict[str, Any]],
    template_execution: Dict[str, Any],
    action_context: Dict[str, Any],
    knowledge_profile=None,
) -> Dict[str, Any]:
    graph_codes = _delta_codes(graph_delta)
    ring_codes = _delta_codes(ring_topology)
    fg_codes = _delta_codes(fg_delta)
    observed = graph_codes | ring_codes | fg_codes
    profile = _knowledge_family_profile(
        cat,
        action_context=action_context,
        knowledge_profile=knowledge_profile,
    )
    declared_family = _family_evidence_summary(action_context.get("family_evidence"), ())
    base_supporting_codes: List[str] = []
    if action_context.get("intended_deltas"):
        base_supporting_codes.append("declared_intended_deltas")
    if action_context.get("mechanistic_evidence"):
        base_supporting_codes.append("declared_mechanistic_evidence")
    if declared_family["provided_keys"]:
        base_supporting_codes.append("declared_family_evidence")
    if template_execution.get("target_in_products"):
        base_supporting_codes.append("forward_template_regenerated_target")

    base_required_evidence: List[str] = []
    if _has_complex_ring_delta(ring_topology) or any(
        code in ring_codes
        for code in ("ring_count_increase", "ring_count_decrease", "ring_bond_edit")
    ):
        base_required_evidence.append("mechanistic topology proof for observed ring changes")
    if action_context.get("preserved_anchors") in (None, "", [], {}):
        base_required_evidence.append("preserved scaffold anchors or handles")
    if action_context.get("changed_bonds") in (None, "", [], {}):
        base_required_evidence.append("declared changed-bond atoms")

    base: Dict[str, Any] = {
        "family_key": "",
        "family_class": "",
        "state": "no_family_interpretation",
        "explained_deltas": [],
        "unexplained_deltas": sorted(observed),
        "forbidden_delta_conflicts": [],
        "supporting_codes": sorted(set(base_supporting_codes)),
        "required_evidence": list(dict.fromkeys(base_required_evidence)),
        "family_evidence": declared_family,
        "policy_effect": {
            "downgrade_codes": [],
            "preserve_codes": list(_PRESERVE_CODES),
        },
        "observed_deltas": {
            "graph": sorted(graph_codes),
            "ring": sorted(ring_codes),
            "fg": sorted(fg_codes),
        },
    }
    if not cat:
        return base
    if profile is None:
        base["state"] = "unregistered_family"
        return base

    allowed = set(profile["allowed_deltas"])
    forbidden = set(profile["forbidden_deltas"])
    explained = observed & allowed
    forbidden_conflicts = observed & forbidden
    family_evidence = _family_evidence_summary(
        action_context.get("family_evidence"),
        profile.get("evidence_keys", []) or [],
    )

    contextual_conflicts: List[str] = []
    if profile["family_class"] == "olefination" and _has_complex_ring_delta(ring_topology):
        if len(precursors) != 1:
            contextual_conflicts.append("intermolecular_olefination_complex_ring_delta")
        elif template_execution and not bool(template_execution.get("target_in_products", False)):
            contextual_conflicts.append("olefination_complex_ring_template_miss")
    if profile["family_class"] == "metathesis" and _has_complex_ring_delta(ring_topology) and len(precursors) != 1:
        contextual_conflicts.append("fragmented_rcm_complex_ring_delta")
    if "apparent_scaffold_jump" in graph_codes:
        contextual_conflicts.append("apparent_scaffold_jump")
    if (
        profile["family_class"] == "small_ring_formation"
        and ring_codes
        & {
            "new_fused_ring_system",
            "new_spiro_center",
            "new_bridged_system",
            "new_medium_ring",
            "new_macrocycle",
        }
    ):
        contextual_conflicts.append("small_ring_family_with_high_risk_topology")

    required_evidence = list(profile["required_evidence"])
    if action_context.get("preserved_anchors") in (None, "", [], {}):
        required_evidence.append("preserved scaffold anchors or handles")
    if action_context.get("changed_bonds") in (None, "", [], {}):
        required_evidence.append("declared changed-bond atoms")

    supporting_codes = [f"family_expected_{profile['family_class']}"]
    if action_context.get("intended_deltas"):
        supporting_codes.append("declared_intended_deltas")
    if action_context.get("mechanistic_evidence"):
        supporting_codes.append("declared_mechanistic_evidence")
    if family_evidence["provided_keys"]:
        supporting_codes.append("declared_family_evidence")
    if template_execution.get("target_in_products"):
        supporting_codes.append("forward_template_regenerated_target")

    unexplained = observed - explained - forbidden_conflicts
    if contextual_conflicts or forbidden_conflicts:
        state = "family_delta_conflict"
    elif explained and unexplained:
        state = "partially_explains_delta_requires_evidence"
    elif explained:
        state = "explains_delta_requires_evidence"
    elif observed:
        state = "does_not_explain_observed_delta"
    else:
        state = "family_no_observed_delta"

    downgrade_codes: List[str] = []
    if (
        explained
        and not contextual_conflicts
        and not forbidden_conflicts
        and profile["family_class"] != "high_risk_topology"
    ):
        downgrade_codes.append("high_risk_topology_requires_independent_evidence")

    return {
        "family_key": profile["family_key"],
        "family_class": profile["family_class"],
        "state": state,
        "explained_deltas": sorted(explained),
        "unexplained_deltas": sorted(unexplained),
        "forbidden_delta_conflicts": sorted(forbidden_conflicts | set(contextual_conflicts)),
        "supporting_codes": sorted(set(supporting_codes)),
        "required_evidence": list(dict.fromkeys(required_evidence)),
        "family_evidence": family_evidence,
        "policy_effect": {
            "downgrade_codes": downgrade_codes,
            "preserve_codes": list(_PRESERVE_CODES),
        },
        "observed_deltas": base["observed_deltas"],
        "risk_level": profile["risk_level"],
        "knowledge_ref": profile["knowledge_ref"],
    }


def validate_reaction_family(
    precursors: Sequence[str],
    target: str,
    reaction_category: Optional[str] = None,
    ring_topology: Optional[Dict[str, Any]] = None,
    graph_delta: Optional[Dict[str, Any]] = None,
    fg_delta: Optional[Dict[str, Any]] = None,
    template_execution: Optional[Dict[str, Any]] = None,
    action_context: Optional[Dict[str, Any]] = None,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """Check reaction-family compatibility with topology deltas.

    ``graph_delta`` and ``fg_delta`` are accepted now so later graph/FG audit
    modules can plug into the same family registry without another public
    signature change.
    """
    del target

    cat = _family_category(reaction_category)
    ring_topology = ring_topology or {}
    template_execution = template_execution or {}
    action_context = action_context or {}
    family_interpretation = _build_family_interpretation(
        cat=cat,
        precursors=precursors,
        ring_topology=ring_topology,
        graph_delta=graph_delta,
        fg_delta=fg_delta,
        template_execution=template_execution,
        action_context=action_context,
        knowledge_profile=knowledge_profile,
    )
    if not cat:
        return {
            "pass": True,
            "has_hard_fail": False,
            "violations": [],
            "summary": "no reaction category, skipped",
            "family_interpretation": family_interpretation,
        }

    delta_codes = set(ring_topology.get("delta_codes", []) or [])
    complex_delta = _has_complex_ring_delta(ring_topology)
    violations: List[Dict[str, Any]] = []

    olefination_keys = {"hwe", "horner", "wadsworth", "emmons", "wittig", "julia", "peterson"}
    is_olefination = any(key in cat for key in olefination_keys)
    is_rcm = "rcm" in cat or "metathesis" in cat or "ring_closing_metathesis" in cat

    if is_olefination and complex_delta:
        if len(precursors) != 1:
            violations.append(make_finding(
                code="olefination_lacks_intramolecular_tether_proof",
                severity="requires_evidence",
                source="reaction_family",
                message=(
                    "olefination family cannot justify a new complex ring system from separate "
                    "intermolecular precursors without intramolecular tether proof"
                ),
                evidence={
                    "precursor_count": len(precursors),
                    "delta_codes": sorted(delta_codes),
                    "template_target_in_products": template_execution.get("target_in_products"),
                    "action_in_ring": action_context.get("in_ring"),
                },
                required_evidence=[
                    "single-precursor intramolecular tether proof",
                    "atom-mapped source of new ring junction atoms",
                ],
            ))
        elif template_execution and not bool(template_execution.get("target_in_products", False)):
            violations.append(make_finding(
                code="olefination_ring_closure_template_miss",
                severity="requires_evidence",
                source="reaction_family",
                message="olefination template did not regenerate the complex ring product",
                evidence={
                    "delta_codes": sorted(delta_codes),
                    "best_match": template_execution.get("best_match", ""),
                    "tanimoto_to_target": template_execution.get("tanimoto_to_target"),
                },
                required_evidence=[
                    "forward template target regeneration",
                    "independent atom-mapped ring-closure proof",
                ],
            ))

    if is_rcm and complex_delta:
        if len(precursors) != 1:
            violations.append(make_finding(
                code="rcm_requires_single_diene_precursor",
                severity="requires_evidence",
                source="reaction_family",
                message="RCM ring construction requires one tethered diene precursor, not disconnected fragments",
                evidence={
                    "precursor_count": len(precursors),
                    "delta_codes": sorted(delta_codes),
                },
                required_evidence=[
                    "single tethered diene precursor",
                    "mapped alkene termini in the same precursor",
                ],
            ))

    # Generic safeguard: a ring-bond template miss on a critical ring delta is
    # not enough by itself to hard-block every family, but it must require proof.
    if (
        action_context.get("in_ring")
        and ring_topology.get("risk_level") == "critical"
        and template_execution
        and not bool(template_execution.get("target_in_products", False))
        and not violations
    ):
        violations.append(make_finding(
            code="high_risk_topology_requires_independent_evidence",
            severity="requires_evidence",
            source="reaction_family",
            message="high-risk topology change requires independent mechanism and atom-source evidence",
            evidence={
                "delta_codes": sorted(delta_codes),
                "best_match": template_execution.get("best_match", ""),
                "tanimoto_to_target": template_execution.get("tanimoto_to_target"),
            },
            required_evidence=[
                "forward template target regeneration",
                "mechanistic topology proof",
            ],
        ))

    hard_fail = any(v.get("severity") == "hard_fail" for v in violations)
    return {
        "pass": not hard_fail,
        "has_hard_fail": hard_fail,
        "violations": violations,
        "summary": "ok" if not violations else f"{len(violations)} ring-family violation(s)",
        "family_interpretation": family_interpretation,
    }
