"""LLM-facing validation evidence packets.

The validation gate decides commit policy. This packet is different: it gives
the LLM a compact, structured explanation of observed graph/FG/topology facts,
declared intent, and missing proof obligations.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping


_INTENT_KEYS = (
    "intended_deltas",
    "expected_ring_change",
    "changed_bonds",
    "preserved_anchors",
    "mechanistic_evidence",
    "family_evidence",
    "rationale_summary",
)


def _dedupe(values: Iterable[Any]) -> List[str]:
    out: List[str] = []
    seen = set()
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def _compact_value(value: Any, *, depth: int = 0) -> Any:
    if depth > 2:
        return str(value)
    if isinstance(value, Mapping):
        compact: Dict[str, Any] = {}
        for idx, (key, item) in enumerate(value.items()):
            if idx >= 8:
                compact["..."] = "truncated"
                break
            compact[str(key)] = _compact_value(item, depth=depth + 1)
        return compact
    if isinstance(value, list):
        return [_compact_value(item, depth=depth + 1) for item in value[:8]]
    if isinstance(value, tuple):
        return [_compact_value(item, depth=depth + 1) for item in list(value)[:8]]
    return value


def _compact_atom_mapping_finding_evidence(evidence: Mapping[str, Any]) -> Dict[str, Any]:
    compact: Dict[str, Any] = {}
    for key in (
        "status",
        "confidence",
        "mcs_product_coverage",
        "mcs_num_atoms",
        "product_heavy_atoms",
        "missing_fields",
        "product_match_count",
        "precursor_match_count",
    ):
        value = evidence.get(key)
        if value not in (None, "", [], {}):
            compact[key] = value
    sample = evidence.get("sample") or []
    if sample:
        compact["sample"] = [
            {
                out_key: row.get(out_key)
                for out_key in (
                    "product_atoms",
                    "product_symbols",
                    "precursor_indices",
                    "precursor_local_atoms",
                    "same_precursor",
                    "product_bond_in_ring",
                )
                if row.get(out_key) not in (None, "", [], {})
            }
            for row in list(sample)[:2]
            if isinstance(row, Mapping)
        ]
    return compact


def _finding_row(finding: Mapping[str, Any]) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "code": str(finding.get("code", "") or ""),
        "severity": str(finding.get("severity", "") or ""),
        "source": str(finding.get("source", "") or ""),
        "message": str(finding.get("message", "") or ""),
    }
    evidence = finding.get("evidence") or {}
    if evidence:
        if row["source"] == "atom_mapping":
            row["evidence"] = _compact_atom_mapping_finding_evidence(evidence)
        else:
            row["evidence"] = _compact_value(evidence)
    required = finding.get("required_evidence") or []
    if required:
        row["required_evidence"] = _dedupe(required)
    return row


def _collect_findings(checks: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for check_name in ("graph_delta", "ring_topology", "atom_mapping", "fg_delta", "ring_family"):
        payload = checks.get(check_name) or {}
        for finding in (payload.get("findings", []) or []) + (payload.get("violations", []) or []):
            if isinstance(finding, Mapping):
                if (
                    check_name == "atom_mapping"
                    and finding.get("severity") not in {"missing_evidence", "requires_evidence", "hard_fail"}
                ):
                    continue
                row = _finding_row(finding)
                if row.get("code"):
                    rows.append(row)

    deduped: List[Dict[str, Any]] = []
    seen = set()
    for row in rows:
        key = (row.get("source"), row.get("code"), row.get("severity"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def _collect_gate_codes(gate: Mapping[str, Any], key: str) -> List[str]:
    return _dedupe(
        item.get("code", "")
        for item in gate.get(key, []) or []
        if isinstance(item, Mapping)
    )


def _compact_family_interpretation(ring_family: Mapping[str, Any]) -> Dict[str, Any]:
    family = ring_family.get("family_interpretation") or {}
    if not isinstance(family, Mapping) or not family:
        return {}

    policy_effect = family.get("policy_effect") or {}
    compact: Dict[str, Any] = {
        "family_key": str(family.get("family_key", "") or ""),
        "family_class": str(family.get("family_class", "") or ""),
        "state": str(family.get("state", "") or ""),
        "explained_deltas": _dedupe(family.get("explained_deltas", []) or []),
        "unexplained_deltas": _dedupe(family.get("unexplained_deltas", []) or []),
        "forbidden_delta_conflicts": _dedupe(family.get("forbidden_delta_conflicts", []) or []),
        "supporting_codes": _dedupe(family.get("supporting_codes", []) or []),
        "required_evidence": _dedupe(family.get("required_evidence", []) or []),
        "family_evidence": {
            "state": str((family.get("family_evidence") or {}).get("state", "") or ""),
            "provided_keys": _dedupe((family.get("family_evidence") or {}).get("provided_keys", []) or []),
            "missing_keys": _dedupe((family.get("family_evidence") or {}).get("missing_keys", []) or []),
            "required_keys": _dedupe((family.get("family_evidence") or {}).get("required_keys", []) or []),
            "missing_evidence_codes": _dedupe((family.get("family_evidence") or {}).get("missing_evidence_codes", []) or []),
            "provided": _compact_value((family.get("family_evidence") or {}).get("provided", {}) or {}),
        },
    }
    if family.get("risk_level"):
        compact["risk_level"] = str(family.get("risk_level", "") or "")
    if isinstance(policy_effect, Mapping) and policy_effect.get("downgrade_codes"):
        compact["policy_effect"] = {
            "downgrade_codes": _dedupe(policy_effect.get("downgrade_codes", []) or []),
        }
    return compact


def _compact_policy_preview(assessment: Mapping[str, Any]) -> Dict[str, Any]:
    preview = assessment.get("policy_preview") or {}
    if not isinstance(preview, Mapping) or not preview:
        return {}

    compact: Dict[str, Any] = {}
    for key in (
        "dry_run_only",
        "policy_source",
        "family_key",
        "family_class",
        "family_state",
        "risk_level",
        "current_gate_state",
        "current_commit_policy",
        "preview_gate_state",
        "preview_commit_policy",
        "would_change_gate",
        "would_change_commit_policy",
    ):
        value = preview.get(key)
        if value not in (None, "", [], {}):
            compact[key] = value

    for key in (
        "applied_downgrade_codes",
        "preserved_gate_codes",
        "preview_hard_block_codes",
        "preview_override_allowed_codes",
        "preview_missing_evidence_codes",
        "preview_soft_warning_codes",
    ):
        value = preview.get(key) or []
        if value:
            compact[key] = _dedupe(value)
    return compact


def _compact_policy_adjustments(gate: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for item in gate.get("policy_adjustments", []) or []:
        if not isinstance(item, Mapping):
            continue
        row: Dict[str, Any] = {}
        for key in (
            "activation_scope",
            "policy_source",
            "family_key",
            "family_class",
            "family_state",
            "risk_level",
            "from_bucket",
            "to_bucket",
            "reason",
        ):
            value = item.get(key)
            if value not in (None, "", [], {}):
                row[key] = value
        codes = item.get("applied_downgrade_codes") or []
        if codes:
            row["applied_downgrade_codes"] = _dedupe(codes)
        if row:
            rows.append(row)
    return rows


def _compact_atom_mapping(atom_mapping: Mapping[str, Any]) -> Dict[str, Any]:
    if not isinstance(atom_mapping, Mapping) or not atom_mapping:
        return {}
    compact: Dict[str, Any] = {}
    for key in (
        "status",
        "confidence",
        "mcs_product_coverage",
    ):
        value = atom_mapping.get(key)
        if value not in (None, "", [], {}):
            compact[key] = value
    ambiguity = atom_mapping.get("ambiguity") or {}
    if ambiguity and ambiguity.get("ambiguous"):
        compact["ambiguity"] = {
            "ambiguous": bool(ambiguity.get("ambiguous", False)),
            "product_match_count": ambiguity.get("product_match_count"),
            "precursor_match_count": ambiguity.get("precursor_match_count"),
        }
    declared = atom_mapping.get("declared_evidence") or {}
    if declared:
        missing_fields = declared.get("missing_fields") or []
        if missing_fields or not (
            declared.get("changed_bonds_provided")
            and declared.get("preserved_anchors_provided")
        ):
            compact["declared_evidence"] = {
                "changed_bonds": bool(declared.get("changed_bonds_provided", False)),
                "preserved_anchors": bool(declared.get("preserved_anchors_provided", False)),
                "family_evidence": bool(declared.get("family_evidence_provided", False)),
                "missing_fields": list(missing_fields),
            }
    compact["change_counts"] = {
        "formed": len(atom_mapping.get("formed_bonds", []) or []),
        "cleaved": len(atom_mapping.get("cleaved_bonds", []) or []),
        "order": len(atom_mapping.get("bond_order_changes", []) or []),
    }
    for key in ("formed_bonds", "cleaved_bonds"):
        value = atom_mapping.get(key) or []
        if value:
            compact[key] = [
                {
                    "p": row.get("product_atoms"),
                    "sym": row.get("product_symbols"),
                    "prec": row.get("precursor_indices"),
                    "same": row.get("same_precursor"),
                    "ring": row.get("product_bond_in_ring"),
                }
                for row in list(value)[:2]
                if isinstance(row, Mapping)
            ]
    for key in ("unmapped_product_atoms",):
        value = atom_mapping.get(key)
        if value not in (None, "", [], {}):
            compact[key] = _compact_value(value)
    return compact


def _gate_interpretation(gate: Mapping[str, Any]) -> Dict[str, str]:
    state = str(gate.get("gate_state", "") or "")
    if state == "override_required":
        return {
            "meaning": "proof_obligation_not_disproof",
            "next_step": "revise_or_override",
        }
    if state == "hard_block":
        return {
            "meaning": "blocking_contradiction",
            "next_step": "choose_another_candidate_or_fix_precursors",
        }
    if state == "missing_evidence":
        return {
            "meaning": "inconclusive_missing_evidence",
            "next_step": "commit_only_after_reasoning_or_add_evidence",
        }
    if state == "warning":
        return {
            "meaning": "soft_warning",
            "next_step": "address_warning_in_reasoning",
        }
    if state == "pass":
        return {
            "meaning": "no_gate_objection",
            "next_step": "normal_chemistry_review",
        }
    return {
        "meaning": "unknown_gate_state",
        "next_step": "inspect_validation_payload",
    }


def build_validation_evidence_packet(forward_validation: Mapping[str, Any]) -> Dict[str, Any]:
    """Build a compact packet for LLM reasoning after sandbox validation."""
    checks = forward_validation.get("checks") or {}
    assessment = forward_validation.get("assessment") or forward_validation
    gate = assessment.get("gate") or forward_validation.get("gate") or {}
    action_context = checks.get("action_context") or {}
    graph_delta = checks.get("graph_delta") or {}
    ring_topology = checks.get("ring_topology") or {}
    atom_mapping_raw = checks.get("atom_mapping") or {}
    atom_mapping = _compact_atom_mapping(atom_mapping_raw)
    fg_delta = checks.get("fg_delta") or {}
    template_exec = checks.get("template_execution") or {}
    atom_balance = checks.get("atom_balance") or {}
    scaffold = checks.get("scaffold_alignment") or {}
    ring_family = checks.get("ring_family") or {}
    precursor_state_raw = checks.get("precursor_state") or {}
    family_interpretation = _compact_family_interpretation(ring_family)
    policy_preview = _compact_policy_preview(assessment)
    policy_adjustments = _compact_policy_adjustments(gate)

    declared_intent = {
        key: _compact_value(action_context.get(key))
        for key in _INTENT_KEYS
        if action_context.get(key) not in (None, "", [], {})
    }
    findings = _collect_findings(checks)
    required_evidence = _dedupe(
        req
        for finding in findings
        for req in finding.get("required_evidence", []) or []
    )
    if family_interpretation:
        required_evidence = _dedupe(required_evidence + list(family_interpretation.get("required_evidence", []) or []))

    review_focus: List[str] = []
    contradiction_codes = _dedupe(
        _collect_gate_codes(gate, "hard_blocks")
        + _collect_gate_codes(gate, "override_allowed")
        + [
            finding.get("code", "")
            for finding in findings
            if finding.get("severity") in {"hard_fail", "requires_evidence"}
        ]
    )
    missing_evidence_codes = _dedupe(
        _collect_gate_codes(gate, "missing_evidence")
        + [
            finding.get("code", "")
            for finding in findings
            if finding.get("severity") == "missing_evidence"
        ]
    )
    supporting_codes = _dedupe(
        finding.get("code", "")
        for finding in findings
        if finding.get("severity") in {"info", "warning"}
    )
    if family_interpretation:
        supporting_codes = _dedupe(supporting_codes + list(family_interpretation.get("supporting_codes", []) or []))
        contradiction_codes = _dedupe(
            contradiction_codes + list(family_interpretation.get("forbidden_delta_conflicts", []) or [])
        )
        family_evidence = family_interpretation.get("family_evidence") or {}
        if family_evidence:
            review_focus.append(
                "Family evidence keys: "
                f"provided={', '.join(family_evidence.get('provided_keys', []) or []) or 'none'}; "
                f"missing={', '.join(family_evidence.get('missing_keys', []) or []) or 'none'}."
            )
            missing_evidence_codes = _dedupe(
                missing_evidence_codes + list(family_evidence.get("missing_evidence_codes", []) or [])
            )
    if policy_adjustments:
        for item in policy_adjustments:
            review_focus.append(
                f"Policy adjustment applied: {item.get('family_key', '') or 'family'} "
                f"{item.get('from_bucket', '') or 'override'} -> "
                f"{item.get('to_bucket', '') or 'missing_evidence'}."
            )
    if policy_preview and policy_preview.get("would_change_gate"):
        review_focus.append(
            f"Policy preview: {policy_preview.get('current_gate_state', '') or 'unknown'} -> "
            f"{policy_preview.get('preview_gate_state', '') or 'unknown'}."
        )
    if atom_mapping:
        review_focus.append(
            "Atom-map audit: "
            f"{atom_mapping.get('status', '') or 'unknown'} / "
            f"{atom_mapping.get('confidence', '') or 'unknown'}; "
            f"formed={len(atom_mapping.get('formed_bonds', []) or [])}, "
            f"cleaved={len(atom_mapping.get('cleaved_bonds', []) or [])}."
        )
        missing_fields = ((atom_mapping.get("declared_evidence") or {}).get("missing_fields") or [])
        if missing_fields:
            review_focus.append(
                "Custom atom-source audit fields missing: " + ", ".join(missing_fields) + "."
            )
    if contradiction_codes:
        review_focus.append("Resolve mechanism/family conflict before commit.")
    if required_evidence:
        review_focus.append("Provide the listed required evidence or choose another action.")
    if missing_evidence_codes:
        review_focus.append("Template absence is inconclusive; do not treat it as positive evidence.")
    if not declared_intent and action_context.get("source") in {"custom_precursors", "llm_proposed"}:
        review_focus.append("Custom action lacks declared intended deltas and preserved anchors.")
    if family_interpretation:
        review_focus.append(
            f"Family interpretation: {family_interpretation.get('family_key', '') or 'unregistered'} / "
            f"{family_interpretation.get('state', '') or 'unknown'}."
        )

    precursor_state = {}
    if precursor_state_raw.get("status") not in (None, "", "closed_shell"):
        precursor_state = {
            key: _compact_value(precursor_state_raw.get(key))
            for key in (
                "status",
                "declared",
                "total_radical_electrons",
                "open_shell_precursors",
                "elemental_metal_reagents",
            )
            if precursor_state_raw.get(key) not in (None, "", [], {})
        }

    return {
        "source": "validation_evidence_packet",
        "gate_state": gate.get("gate_state", ""),
        "commit_policy": gate.get("commit_policy", ""),
        "gate_interpretation": _gate_interpretation(gate),
        "reaction_context": {
            "source": action_context.get("source", ""),
            "action_id": action_context.get("action_id", ""),
            "reaction_type": action_context.get("reaction_type", ""),
            "risk_tags": list(action_context.get("risk_tags", []) or []),
        },
        "declared_intent": declared_intent,
        "family_interpretation": family_interpretation,
        "policy_preview": policy_preview,
        "policy_adjustments": policy_adjustments,
        "precursor_state": precursor_state,
        "atom_mapping": atom_mapping,
        "observed_deltas": {
            "graph": list(graph_delta.get("delta_codes", []) or []),
            "ring": list(ring_topology.get("delta_codes", []) or []),
            "atom_mapping": list(atom_mapping_raw.get("delta_codes", []) or []),
            "fg": list(fg_delta.get("delta_codes", []) or []),
        },
        "validation_signals": {
            "template_attempted": bool(template_exec.get("attempted", False)),
            "template_target_in_products": template_exec.get("target_in_products"),
            "template_best_match": template_exec.get("best_match", ""),
            "template_tanimoto_to_target": template_exec.get("tanimoto_to_target"),
            "atom_balance_note": atom_balance.get("note", ""),
            "scaffold_coverage_ratio": scaffold.get("coverage_ratio"),
            "scaffold_summary": scaffold.get("summary", ""),
        },
        "contradiction_codes": contradiction_codes,
        "missing_evidence_codes": missing_evidence_codes,
        "supporting_codes": supporting_codes,
        "required_evidence": required_evidence,
        "key_findings": findings[:12],
        "llm_review_focus": review_focus[:7],
    }
