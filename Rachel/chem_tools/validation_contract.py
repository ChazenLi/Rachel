"""Canonical public validation vocabulary.

The chemistry validators retain their detailed legacy payloads for diagnostics
and historical session loading. Public Rachel commands project those payloads
through this module so execution state, chemistry findings, proof obligations,
and tool limitations are not conflated.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple


SCHEMA_VERSION = "rachel.validation.v2"

_STATE_MAP = {
    "hard_block": "blocked",
    "override_required": "proof_required",
    "missing_evidence": "inconclusive",
    "warning": "warning",
    "pass": "clear",
}

_SYSTEM_ERROR_CODES = {
    "invalid_validation_input",
    "validation_payload_missing",
    "validation_unavailable",
    "validation_gate_missing",
}

_TOOL_LIMIT_CODES = {
    "template_not_attempted",
    "template_target_not_generated",
}

_PUBLIC_TERMS: Dict[str, Tuple[str, str]] = {
    "template_not_attempted": (
        "forward_template_unavailable",
        "No executable forward template was available; this is a tool-coverage limit, not chemical disproof.",
    ),
    "template_target_not_generated": (
        "forward_template_target_not_regenerated",
        "The available forward template did not regenerate the target; independent chemistry evidence is still required.",
    ),
    "atom_mapping_declared_evidence_missing": (
        "declared_atom_source_evidence_missing",
        "The custom action is missing declared changed-bond, atom-source, or preserved-anchor evidence.",
    ),
    "unsupported_new_fused_medium_ring": (
        "new_fused_medium_ring_requires_evidence",
        "A new fused medium-ring topology was detected and requires explicit bond-source and mechanism evidence.",
    ),
    "unsupported_new_fused_ring_system": (
        "new_fused_ring_requires_evidence",
        "A new fused-ring topology was detected and requires explicit bond-source and mechanism evidence.",
    ),
    "unsupported_new_spiro_center": (
        "new_spiro_center_requires_evidence",
        "A new spiro junction was detected and requires explicit junction-source and mechanism evidence.",
    ),
    "unsupported_new_bridged_system": (
        "new_bridged_system_requires_evidence",
        "A new bridged-ring topology was detected and requires explicit bridgehead-source and mechanism evidence.",
    ),
    "unsupported_new_macrocycle": (
        "new_macrocycle_requires_evidence",
        "A new macrocycle was detected and requires explicit tether, atom-source, and cyclization evidence.",
    ),
    "critical_ring_delta_template_miss": (
        "high_risk_topology_requires_independent_evidence",
        "A high-risk topology change needs independent mechanism and atom-source evidence; template absence is not disproof.",
    ),
    "fusion_bond_infeasible": (
        "fusion_bond_mapping_requires_evidence",
        "The mapped fusion bond is not reconciled with the declared event; revise the precursor map or provide explicit bond-source evidence.",
    ),
    "scaffold_not_aligned": (
        "major_scaffold_not_inherited",
        "No single precursor inherits enough of the target scaffold; provide site-specific or scaffold-assembly evidence.",
    ),
    "forbidden_fg": (
        "functional_group_condition_conflict",
        "The proposed conditions may conflict with retained functional groups; resolve conditions, protection, or precursor design.",
    ),
}


def canonical_finding_code(code: Any) -> str:
    legacy_code = str(code or "")
    return _PUBLIC_TERMS.get(legacy_code, (legacy_code, ""))[0]


def canonical_gate_state(state: Any) -> str:
    legacy_state = str(state or "")
    return _STATE_MAP.get(legacy_state, legacy_state or "unknown")


def _dedupe_rows(rows: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    result: List[Dict[str, Any]] = []
    seen = set()
    for row in rows:
        key = (row.get("code"), row.get("source"))
        if key in seen:
            continue
        seen.add(key)
        result.append(row)
    return result


def _gate_from_forward_validation(forward_validation: Mapping[str, Any]) -> Mapping[str, Any]:
    assessment = forward_validation.get("assessment") or {}
    gate = assessment.get("gate") if isinstance(assessment, Mapping) else None
    if not gate:
        gate = forward_validation.get("gate") or forward_validation.get("validation_gate") or {}
    return gate if isinstance(gate, Mapping) else {}


def _public_finding(item: Mapping[str, Any]) -> Dict[str, Any]:
    legacy_code = str(item.get("code", "") or "")
    public_code, public_message = _PUBLIC_TERMS.get(
        legacy_code,
        (legacy_code, str(item.get("message", "") or "")),
    )
    row: Dict[str, Any] = {
        "code": public_code,
        "source": str(item.get("source", "") or ""),
        "message": public_message,
    }
    required = item.get("required_evidence") or []
    if required:
        row["required_evidence"] = [str(value) for value in required if str(value)]
    evidence = item.get("evidence") or {}
    if isinstance(evidence, Mapping) and evidence:
        row["evidence"] = dict(evidence)
    return {key: value for key, value in row.items() if value not in ("", [], {})}


def _gate_rows(gate: Mapping[str, Any], key: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for item in gate.get(key, []) or []:
        if isinstance(item, Mapping) and item.get("code"):
            rows.append(_public_finding(item))
        elif item:
            rows.append(_public_finding({"code": str(item)}))
    return rows


def _template_observation(
    forward_validation: Mapping[str, Any],
    validation_micro: Mapping[str, Any],
) -> Dict[str, Any]:
    checks = forward_validation.get("checks") or {}
    template = checks.get("template_execution") if isinstance(checks, Mapping) else {}
    template = template if isinstance(template, Mapping) else {}
    attempted = bool(
        template.get("attempted", validation_micro.get("template_attempted", False))
    )
    target_in_products = template.get("target_in_products")
    if not attempted:
        return {"status": "not_attempted"}
    status = "target_regenerated" if target_in_products is True else "target_not_regenerated"
    row: Dict[str, Any] = {"status": status}
    similarity = template.get("tanimoto_to_target")
    if similarity is None:
        similarity = validation_micro.get("template_product_similarity")
    if similarity is None:
        similarity = validation_micro.get("template_match")
    if similarity is not None:
        row["product_similarity"] = similarity
    return row


def _mapping_observation(validation_micro: Mapping[str, Any]) -> Dict[str, Any]:
    status = str(validation_micro.get("atom_mapping_status", "") or "")
    if not status:
        return {}
    row: Dict[str, Any] = {
        "method": "mcs_heuristic",
        "status": status,
    }
    coverage = validation_micro.get("atom_mapping_mcs_product_coverage")
    if coverage is not None:
        row["product_coverage"] = coverage
    if validation_micro.get("atom_mapping_ambiguous"):
        row["ambiguous"] = True
    formed = validation_micro.get("atom_mapping_formed_bond_count")
    cleaved = validation_micro.get("atom_mapping_cleaved_bond_count")
    if formed is not None or cleaved is not None:
        row["change_counts"] = {
            "formed": int(formed or 0),
            "cleaved": int(cleaved or 0),
        }
    return row


def _site_observation(site_audit: Mapping[str, Any]) -> Dict[str, Any]:
    if not site_audit:
        return {}
    site_retentive = bool(site_audit.get("site_retentive"))
    if not site_retentive:
        return {"state": "not_applicable"}
    row: Dict[str, Any] = {
        "state": "mapped" if site_audit.get("ok", True) else "unavailable",
        "mode": "strict" if site_audit.get("strict_audit_required") else "informational",
        "same_core": True,
    }
    changed_count = site_audit.get("changed_site_count")
    if changed_count is not None:
        row["changed_site_count"] = changed_count
    changed_sites = site_audit.get("changed_sites") or []
    if changed_sites:
        row["changed_sites"] = list(changed_sites)[:3]
    summary = str(site_audit.get("summary", "") or "")
    if summary:
        row["summary"] = summary
    return row


def _eas_site_observation(eas_site_audit: Mapping[str, Any]) -> Dict[str, Any]:
    if not eas_site_audit or not eas_site_audit.get("applicable"):
        return {}
    directors = [
        {
            "rule_id": str(item.get("rule_id", "") or ""),
            "ring_atom_idx": item.get("ring_atom_idx"),
            "directing_positions": list(item.get("directing_positions", []) or []),
            "electronic_effect": str(item.get("electronic_effect", "") or ""),
            "knowledge_ref": dict(item.get("knowledge_ref", {}) or {}),
        }
        for item in eas_site_audit.get("directors", []) or []
        if isinstance(item, Mapping) and item.get("rule_id")
    ]
    row: Dict[str, Any] = {
        "state": str(eas_site_audit.get("state", "") or ""),
        "method": str(eas_site_audit.get("method", "") or ""),
        "reason_code": str(eas_site_audit.get("reason_code", "") or ""),
        "atom_mapping_is_selectivity_proof": False,
        "desired_atom_idx": eas_site_audit.get("desired_atom_idx"),
        "directors": directors,
        "candidate_sites": list(eas_site_audit.get("candidate_sites", []) or []),
        "summary": str(eas_site_audit.get("summary", "") or ""),
    }
    return {key: value for key, value in row.items() if value not in ("", [], {})}


def _legacy_evidence_packet(
    forward_validation: Mapping[str, Any],
    validation_micro: Mapping[str, Any],
) -> Mapping[str, Any]:
    assessment = forward_validation.get("assessment") or {}
    candidates = [
        assessment.get("evidence_packet") if isinstance(assessment, Mapping) else None,
        forward_validation.get("evidence_packet"),
        validation_micro.get("evidence_packet"),
    ]
    for candidate in candidates:
        if isinstance(candidate, Mapping) and candidate:
            return candidate
    return {}


def _mechanism_interpretation(packet: Mapping[str, Any]) -> Dict[str, Any]:
    family = packet.get("family_interpretation") or {}
    if not isinstance(family, Mapping) or not family:
        return {}
    label = str(family.get("family_key", "") or "")
    state = str(family.get("state", "") or "")
    family_evidence = family.get("family_evidence") or {}
    if state in {"unregistered_family", "no_family_interpretation"}:
        return {}
    if not label and not state and not family_evidence:
        return {}
    row: Dict[str, Any] = {
        "label": label,
        "class": str(family.get("family_class", "") or ""),
        "state": state,
        "explained_observations": list(family.get("explained_deltas", []) or []),
        "unresolved_observations": list(family.get("unexplained_deltas", []) or []),
    }
    if isinstance(family_evidence, Mapping) and family_evidence:
        row["evidence"] = {
            "state": str(family_evidence.get("state", "") or ""),
            "provided_keys": list(family_evidence.get("provided_keys", []) or []),
            "missing_keys": list(family_evidence.get("missing_keys", []) or []),
        }
    return {key: value for key, value in row.items() if value not in ("", [], {})}


def _policy_adjustments(packet: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for item in packet.get("policy_adjustments", []) or []:
        if not isinstance(item, Mapping):
            continue
        codes = [
            canonical_finding_code(code)
            for code in item.get("applied_downgrade_codes", []) or []
            if str(code)
        ]
        row = {
            "basis": "mechanism_interpretation",
            "from_category": str(item.get("from_bucket", "") or ""),
            "to_category": str(item.get("to_bucket", "") or ""),
            "finding_codes": codes,
        }
        rows.append({key: value for key, value in row.items() if value not in ("", [], {})})
    return rows


def build_validation_contract(
    forward_validation: Optional[Mapping[str, Any]],
    *,
    validation_micro: Optional[Mapping[str, Any]] = None,
    site_audit: Optional[Mapping[str, Any]] = None,
    eas_site_audit: Optional[Mapping[str, Any]] = None,
    execution_success: Optional[bool] = None,
    execution_scope: str = "sandbox_validation",
) -> Dict[str, Any]:
    """Project legacy validation payloads into the canonical public contract."""
    forward_validation = forward_validation if isinstance(forward_validation, Mapping) else {}
    validation_micro = validation_micro if isinstance(validation_micro, Mapping) else {}
    site_audit = site_audit if isinstance(site_audit, Mapping) else {}
    eas_site_audit = (
        eas_site_audit if isinstance(eas_site_audit, Mapping) else {}
    )
    gate = _gate_from_forward_validation(forward_validation)
    legacy_packet = _legacy_evidence_packet(forward_validation, validation_micro)

    hard_rows = _gate_rows(gate, "hard_blocks")
    system_errors = [
        row for row in hard_rows if row.get("code") in _SYSTEM_ERROR_CODES
    ]
    contradictions = [
        row for row in hard_rows if row.get("code") not in _SYSTEM_ERROR_CODES
    ]
    proof_obligations = _gate_rows(gate, "override_allowed")
    missing_rows = _gate_rows(gate, "missing_evidence")
    tool_limits = [
        row
        for row in missing_rows
        if row.get("code") in {
            _PUBLIC_TERMS[code][0] for code in _TOOL_LIMIT_CODES
        }
    ]
    evidence_gaps = [row for row in missing_rows if row not in tool_limits]
    warnings = _gate_rows(gate, "soft_warnings")

    legacy_state = str(gate.get("gate_state", gate.get("state", "")) or "")
    decision_gate: Dict[str, Any] = {
        "state": canonical_gate_state(legacy_state),
    }
    if decision_gate["state"] == "blocked":
        decision_gate["block_type"] = (
            "system_error" if system_errors and not contradictions else "chemical_contradiction"
        )

    observations: Dict[str, Any] = {
        "template": _template_observation(forward_validation, validation_micro),
    }
    mapping = _mapping_observation(validation_micro)
    if not mapping and isinstance(legacy_packet.get("atom_mapping"), Mapping):
        legacy_mapping = legacy_packet.get("atom_mapping") or {}
        mapping = {
            "method": "mcs_heuristic",
            "status": str(legacy_mapping.get("status", "") or ""),
        }
        if legacy_mapping.get("mcs_product_coverage") is not None:
            mapping["product_coverage"] = legacy_mapping.get("mcs_product_coverage")
        if legacy_mapping.get("ambiguity"):
            mapping["ambiguous"] = True
        if legacy_mapping.get("change_counts"):
            mapping["change_counts"] = dict(legacy_mapping.get("change_counts") or {})
        mapping = {key: value for key, value in mapping.items() if value not in ("", [], {})}
    if mapping:
        observations["atom_mapping"] = mapping
    graph_codes = list(validation_micro.get("graph_delta_codes", []) or [])
    ring_codes = list(validation_micro.get("ring_delta_codes", []) or [])
    fg_codes = list(validation_micro.get("fg_delta_codes", []) or [])
    if graph_codes:
        observations["graph_deltas"] = graph_codes
    if ring_codes:
        observations["ring_deltas"] = ring_codes
    if fg_codes:
        observations["functional_group_deltas"] = fg_codes
    packet_deltas = legacy_packet.get("observed_deltas") or {}
    if isinstance(packet_deltas, Mapping):
        if not graph_codes and packet_deltas.get("graph"):
            observations["graph_deltas"] = list(packet_deltas.get("graph") or [])
        if not ring_codes and packet_deltas.get("ring"):
            observations["ring_deltas"] = list(packet_deltas.get("ring") or [])
        if not fg_codes and packet_deltas.get("fg"):
            observations["functional_group_deltas"] = list(packet_deltas.get("fg") or [])
    site = _site_observation(site_audit)
    if site:
        observations["site_fidelity"] = site
    eas_site = _eas_site_observation(eas_site_audit)
    if eas_site:
        observations["eas_site_selectivity"] = eas_site
    precursor_state = (
        (forward_validation.get("checks") or {}).get("precursor_state")
        or legacy_packet.get("precursor_state")
        or {}
    )
    if precursor_state.get("status") == "open_shell_detected":
        observations["precursor_state"] = {
            "status": "open_shell_detected",
            "declared": bool(precursor_state.get("declared")),
            "total_radical_electrons": precursor_state.get(
                "total_radical_electrons", 0
            ),
            "open_shell_precursors": precursor_state.get(
                "open_shell_precursors", []
            ),
        }
        if precursor_state.get("elemental_metal_reagents"):
            observations["precursor_state"]["elemental_metal_reagents"] = precursor_state.get(
                "elemental_metal_reagents", []
            )
    elif precursor_state.get("status") == "elemental_metal_reagent":
        observations["precursor_state"] = {
            "status": "elemental_metal_reagent",
            "elemental_metal_reagents": precursor_state.get(
                "elemental_metal_reagents", []
            ),
        }

    declared_action = legacy_packet.get("declared_intent") or {}
    if not isinstance(declared_action, Mapping):
        declared_action = {}
    mechanism = _mechanism_interpretation(legacy_packet)
    adjustments = _policy_adjustments(legacy_packet)
    if adjustments:
        decision_gate["policy_adjustments"] = adjustments

    next_step = {
        "blocked": "revise_candidate_or_rerun_validator",
        "proof_required": "add_evidence_revise_or_override",
        "inconclusive": "complete_chemistry_review",
        "warning": "address_warnings",
        "clear": "normal_chemistry_review",
    }.get(decision_gate["state"], "inspect_validation")

    contract: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "execution": {
            "status": "completed" if execution_success else "failed"
            if execution_success is not None
            else "unknown",
            "scope": execution_scope,
        },
        "decision_gate": {
            key: value for key, value in decision_gate.items() if value not in ("", [], {})
        },
        "contradictions": _dedupe_rows(contradictions),
        "proof_obligations": _dedupe_rows(proof_obligations),
        "evidence_gaps": _dedupe_rows(evidence_gaps),
        "tool_limits": _dedupe_rows(tool_limits),
        "warnings": _dedupe_rows(warnings),
        "system_errors": _dedupe_rows(system_errors),
        "observations": observations,
        "declared_action": dict(declared_action),
        "mechanism_interpretation": mechanism,
        "recommended_next_step": next_step,
    }
    return {
        key: value
        for key, value in contract.items()
        if value not in ("", [], {}) or key in {"contradictions", "proof_obligations"}
    }
