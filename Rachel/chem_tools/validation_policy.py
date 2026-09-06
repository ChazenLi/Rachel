"""Validation gate policy aggregation.

Audit modules report facts. This module is the commit-facing policy layer that
aggregates atom balance, template execution, topology, family, and FG evidence
into hard_block / override_required / missing_evidence / warning / pass.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional


def _validation_gate_item(
    code: str,
    source: str,
    message: str,
    evidence: Optional[Dict[str, Any]] = None,
    required_evidence: Optional[List[str]] = None,
) -> Dict[str, Any]:
    item: Dict[str, Any] = {
        "code": code,
        "source": source,
        "message": message,
    }
    if evidence:
        item["evidence"] = evidence
    if required_evidence:
        item["required_evidence"] = required_evidence
    return item


def _missing_protecting_group_source(checks: Mapping[str, Any]) -> Dict[str, Any]:
    atom_balance = checks.get("atom_balance") or {}
    fg_delta = checks.get("fg_delta") or {}
    graph_delta = checks.get("graph_delta") or {}
    protecting_group_delta = fg_delta.get("protecting_group_delta") or {}
    increased = protecting_group_delta.get("increased") or {}
    if not increased:
        return {}
    if "protecting_group_added" not in set(fg_delta.get("delta_codes", []) or []):
        return {}
    if list(graph_delta.get("delta_codes", []) or []):
        return {}
    if float(graph_delta.get("mcs_precursor_coverage", 0.0) or 0.0) < 0.99:
        return {}
    if not (
        atom_balance.get("skeleton_imbalance")
        or atom_balance.get("severe_imbalance")
    ):
        return {}
    return {
        "protecting_groups": dict(increased),
        "missing_atoms": dict(atom_balance.get("adjusted_excess", {}) or {}),
        "atom_balance_note": atom_balance.get("note", ""),
    }


def _collect_gate_items(gate: Mapping[str, Any], key: str) -> List[Dict[str, Any]]:
    return [
        item
        for item in gate.get(key, []) or []
        if isinstance(item, Mapping)
    ]


def _collect_gate_codes(gate: Mapping[str, Any], key: str) -> List[str]:
    seen = set()
    codes: List[str] = []
    for item in _collect_gate_items(gate, key):
        code = str(item.get("code", "") or "").strip()
        if code and code not in seen:
            seen.add(code)
            codes.append(code)
    return codes


_STAGE4_DOWNGRADE_CODE = "high_risk_topology_requires_independent_evidence"
_STAGE4_CONTROLLED_FAMILY_KEYS = {
    "small_ring_formation",
    "ring_opening",
    "intramolecular_ring_closure",
}
_STAGE4_HIGH_RISK_RING_CODES = {
    "new_fused_ring_system",
    "new_spiro_center",
    "new_bridged_system",
    "new_medium_ring",
    "new_macrocycle",
}


def _family_interpretation(forward_validation: Mapping[str, Any]) -> Mapping[str, Any]:
    checks = forward_validation.get("checks") or {}
    ring_family = checks.get("ring_family") or {}
    family = ring_family.get("family_interpretation") or {}
    return family if isinstance(family, Mapping) else {}


def _family_observed_codes(family: Mapping[str, Any], bucket: str) -> set[str]:
    observed = family.get("observed_deltas") or {}
    if not isinstance(observed, Mapping):
        return set()
    return {str(code) for code in (observed.get(bucket, []) or []) if str(code)}


def _stage4_controlled_scope(
    forward_validation: Mapping[str, Any],
    family: Mapping[str, Any],
) -> str:
    family_key = str(family.get("family_key", "") or "")
    state = str(family.get("state", "") or "")
    if family_key not in _STAGE4_CONTROLLED_FAMILY_KEYS:
        return ""
    if family.get("forbidden_delta_conflicts"):
        return ""

    ring_codes = _family_observed_codes(family, "ring")
    if ring_codes & _STAGE4_HIGH_RISK_RING_CODES:
        return ""

    if family_key == "small_ring_formation":
        if state != "explains_delta_requires_evidence":
            return ""
        if "new_strained_small_ring" not in ring_codes:
            return ""
        return "stage4_controlled_family_downgrade"

    if family_key == "ring_opening":
        if state not in {
            "explains_delta_requires_evidence",
            "partially_explains_delta_requires_evidence",
        }:
            return ""
        if "ring_count_decrease" not in ring_codes or "ring_count_increase" in ring_codes:
            return ""
        return "stage4_controlled_family_downgrade"

    if family_key == "intramolecular_ring_closure":
        if state != "explains_delta_requires_evidence":
            return ""
        precursors = forward_validation.get("precursors") or []
        if isinstance(precursors, list) and len(precursors) != 1:
            return ""
        if "ring_count_increase" not in ring_codes:
            return ""
        return "stage4_controlled_family_downgrade"

    return ""


def _apply_controlled_family_downgrades(
    forward_validation: Mapping[str, Any],
    *,
    hard_blocks: List[Dict[str, Any]],
    override_allowed: List[Dict[str, Any]],
    missing_evidence: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    if hard_blocks or not override_allowed:
        return []

    family = _family_interpretation(forward_validation)
    activation_scope = _stage4_controlled_scope(forward_validation, family)
    if not activation_scope:
        return []

    policy_effect = family.get("policy_effect") or {}
    if not isinstance(policy_effect, Mapping):
        return []
    downgrade_codes = {
        str(code)
        for code in (policy_effect.get("downgrade_codes", []) or [])
        if str(code)
    }
    preserve_codes = {
        str(code)
        for code in (policy_effect.get("preserve_codes", []) or [])
        if str(code)
    }
    if _STAGE4_DOWNGRADE_CODE not in downgrade_codes:
        return []

    applied_items: List[Dict[str, Any]] = []
    retained_override: List[Dict[str, Any]] = []
    for item in override_allowed:
        code = str(item.get("code", "") or "")
        if code == _STAGE4_DOWNGRADE_CODE and code not in preserve_codes:
            applied_items.append(dict(item))
        else:
            retained_override.append(item)
    if not applied_items:
        return []

    override_allowed[:] = retained_override
    missing_evidence.extend(dict(item) for item in applied_items)

    return [
        {
            "activation_scope": activation_scope,
            "policy_source": "family_interpretation",
            "family_key": str(family.get("family_key", "") or ""),
            "family_class": str(family.get("family_class", "") or ""),
            "family_state": str(family.get("state", "") or ""),
            "risk_level": str(family.get("risk_level", "") or ""),
            "from_bucket": "override_allowed",
            "to_bucket": "missing_evidence",
            "applied_downgrade_codes": _collect_gate_codes(
                {"override_allowed": applied_items}, "override_allowed"
            ),
            "reason": (
                "family interpretation explains this lower-risk topology delta; "
                "keep the finding as missing evidence for LLM review"
            ),
        }
    ]


def _build_policy_preview(
    forward_validation: Mapping[str, Any],
    gate: Mapping[str, Any],
) -> Dict[str, Any]:
    checks = forward_validation.get("checks") or {}
    ring_family = checks.get("ring_family") or {}
    family = ring_family.get("family_interpretation") or {}
    policy_effect = family.get("policy_effect") or {}
    downgrade_codes = {
        str(code)
        for code in (policy_effect.get("downgrade_codes", []) or [])
        if str(code)
    }
    preserve_codes = {
        str(code)
        for code in (policy_effect.get("preserve_codes", []) or [])
        if str(code)
    }

    current_hard_blocks = _collect_gate_items(gate, "hard_blocks")
    current_override_allowed = _collect_gate_items(gate, "override_allowed")
    current_missing_evidence = _collect_gate_items(gate, "missing_evidence")
    current_soft_warnings = _collect_gate_items(gate, "soft_warnings")

    applied_downgrade_items: List[Dict[str, Any]] = []
    preview_override_allowed: List[Dict[str, Any]] = []
    for item in current_override_allowed:
        code = str(item.get("code", "") or "")
        if code and code in downgrade_codes and code not in preserve_codes:
            applied_downgrade_items.append(dict(item))
            continue
        preview_override_allowed.append(dict(item))

    downgraded_bucket = "missing_evidence"
    if not current_missing_evidence:
        downgraded_bucket = "warning"

    preview_missing_evidence = [dict(item) for item in current_missing_evidence]
    preview_soft_warnings = [dict(item) for item in current_soft_warnings]
    if downgraded_bucket == "missing_evidence":
        preview_missing_evidence.extend(dict(item) for item in applied_downgrade_items)
    else:
        preview_soft_warnings.extend(dict(item) for item in applied_downgrade_items)

    if current_hard_blocks:
        preview_gate_state = "hard_block"
        preview_commit_policy = "block"
    elif preview_override_allowed:
        preview_gate_state = "override_required"
        preview_commit_policy = "requires_override"
    elif preview_missing_evidence:
        preview_gate_state = "missing_evidence"
        preview_commit_policy = "allow_with_note"
    elif preview_soft_warnings:
        preview_gate_state = "warning"
        preview_commit_policy = "allow_with_warning"
    else:
        preview_gate_state = "pass"
        preview_commit_policy = "allow"

    current_hard_block_codes = _collect_gate_codes(gate, "hard_blocks")
    current_override_allowed_codes = _collect_gate_codes(gate, "override_allowed")
    current_missing_evidence_codes = _collect_gate_codes(gate, "missing_evidence")
    current_soft_warning_codes = _collect_gate_codes(gate, "soft_warnings")
    preview_hard_block_codes = current_hard_block_codes
    preview_override_allowed_codes = _collect_gate_codes(
        {"override_allowed": preview_override_allowed}, "override_allowed"
    )
    preview_missing_evidence_codes = _collect_gate_codes(
        {"missing_evidence": preview_missing_evidence}, "missing_evidence"
    )
    preview_soft_warning_codes = _collect_gate_codes(
        {"soft_warnings": preview_soft_warnings}, "soft_warnings"
    )
    preserved_gate_codes = sorted(
        code
        for code in (
            current_hard_block_codes
            + current_override_allowed_codes
            + current_missing_evidence_codes
            + current_soft_warning_codes
        )
        if code in preserve_codes
    )

    preview = {
        "dry_run_only": True,
        "policy_source": "family_interpretation",
        "family_key": str(family.get("family_key", "") or ""),
        "family_class": str(family.get("family_class", "") or ""),
        "family_state": str(family.get("state", "") or ""),
        "risk_level": str(family.get("risk_level", "") or ""),
        "current_gate_state": str(gate.get("gate_state", "") or ""),
        "current_commit_policy": str(gate.get("commit_policy", "") or ""),
        "preview_gate_state": preview_gate_state,
        "preview_commit_policy": preview_commit_policy,
        "applied_downgrade_codes": _collect_gate_codes(
            {"override_allowed": applied_downgrade_items}, "override_allowed"
        ),
        "preserved_gate_codes": preserved_gate_codes,
        "preview_hard_block_codes": preview_hard_block_codes,
        "preview_override_allowed_codes": preview_override_allowed_codes,
        "preview_missing_evidence_codes": preview_missing_evidence_codes,
        "preview_soft_warning_codes": preview_soft_warning_codes,
        "would_change_gate": (
            preview_gate_state != str(gate.get("gate_state", "") or "")
            or current_hard_block_codes != preview_hard_block_codes
            or current_override_allowed_codes != preview_override_allowed_codes
            or current_missing_evidence_codes != preview_missing_evidence_codes
            or current_soft_warning_codes != preview_soft_warning_codes
        ),
    }
    preview["would_change_commit_policy"] = (
        preview_commit_policy != str(gate.get("commit_policy", "") or "")
    )
    if applied_downgrade_items:
        preview["applied_downgrade_items"] = applied_downgrade_items
    return preview


def build_validation_policy_preview(forward_validation: Dict[str, Any]) -> Dict[str, Any]:
    """Build a dry-run policy projection without changing the commit gate."""
    if not isinstance(forward_validation, dict):
        return {}

    assessment = forward_validation.get("assessment") or {}
    gate = assessment.get("gate") or forward_validation.get("gate") or {}
    if not isinstance(gate, dict):
        return {}

    preview = _build_policy_preview(forward_validation, gate)
    if not preview.get("family_key"):
        return {}
    return preview


def build_validation_gate(
    forward_validation: Dict[str, Any],
    site_audit: Optional[Dict[str, Any]] = None,
    ring_topology_audit: Optional[Dict[str, Any]] = None,
    action_context: Optional[Dict[str, Any]] = None,
    eas_site_audit: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Classify validation evidence into commit-facing gate categories.

    This is intentionally conservative. Template execution misses and missing
    template IDs are inconclusive evidence, not chemical hard blocks.
    """
    hard_blocks: List[Dict[str, Any]] = []
    soft_warnings: List[Dict[str, Any]] = []
    override_allowed: List[Dict[str, Any]] = []
    missing_evidence: List[Dict[str, Any]] = []

    if not isinstance(forward_validation, dict):
        return {
            "gate_state": "missing_evidence",
            "commit_policy": "allow_with_note",
            "hard_blocks": [],
            "soft_warnings": [],
            "override_allowed": [],
            "missing_evidence": [
                _validation_gate_item(
                    "validation_payload_missing",
                    "validation",
                    "forward_validation payload is missing or not a dict",
                )
            ],
            "llm_override_allowed": False,
            "recommended_action": "Treat validation as inconclusive; gather sandbox evidence before commit.",
        }

    if forward_validation.get("ok") is False:
        hard_blocks.append(
            _validation_gate_item(
                "invalid_validation_input",
                "validation",
                str(forward_validation.get("error") or "validation input failed"),
                {"input": forward_validation.get("input", "")},
            )
        )

    checks = forward_validation.get("checks") or {}
    assessment = forward_validation.get("assessment") or forward_validation
    ring_topology = ring_topology_audit or checks.get("ring_topology") or {}
    action_context = action_context or checks.get("action_context") or {}
    eas_site = eas_site_audit or checks.get("eas_site_selectivity") or {}

    atom_bal = checks.get("atom_balance") or {}
    protecting_group_source = _missing_protecting_group_source(checks)
    if protecting_group_source:
        override_allowed.append(
            _validation_gate_item(
                "protecting_group_source_required",
                "atom_balance",
                "a protecting group appears in the product but its reagent source is absent from the precursor set; add the complete protecting-group donor and revalidate",
                protecting_group_source,
                required_evidence=[
                    "closed-shell protecting-group donor such as the matching silyl, acyl, benzyl, or sulfonyl reagent",
                    "atom balance after adding that donor",
                    "site-selectivity and compatibility rationale for the protected functional group",
                ],
            )
        )
    if atom_bal.get("skeleton_imbalance") and not protecting_group_source:
        hard_blocks.append(
            _validation_gate_item(
                "skeleton_imbalance",
                "atom_balance",
                "product has more skeleton atoms (C/N/S) than the precursor set",
                {"note": atom_bal.get("note", "")},
            )
        )
    if atom_bal.get("severe_imbalance") and not protecting_group_source:
        hard_blocks.append(
            _validation_gate_item(
                "severe_imbalance",
                "atom_balance",
                "product-side non-H atom excess is too large for a one-step transform",
                {"note": atom_bal.get("note", "")},
            )
        )
    if (
        atom_bal
        and atom_bal.get("balanced") is False
        and not atom_bal.get("skeleton_imbalance")
        and not atom_bal.get("severe_imbalance")
    ):
        reagents = list(action_context.get("reagents", []) or [])
        reagent_only_remainder = bool(reagents) and not atom_bal.get("adjusted_excess")
        soft_warnings.append(
            _validation_gate_item(
                "reagent_byproducts_unmodeled" if reagent_only_remainder else "atom_balance_unresolved",
                "atom_balance",
                (
                    "declared current-step reagents leave unmodeled byproducts; product atom sources are present"
                    if reagent_only_remainder
                    else "atom balance has unexplained remainder, but not a hard skeleton/severe block"
                ),
                {
                    "note": atom_bal.get("note", ""),
                    "reagents": reagents,
                },
            )
        )

    template_exec = checks.get("template_execution") or {}
    if template_exec:
        if not bool(template_exec.get("attempted", False)):
            missing_evidence.append(
                _validation_gate_item(
                    "template_not_attempted",
                    "template_execution",
                    "no executable forward template was available; this is missing evidence, not a hard block",
                    {
                        "template_id": template_exec.get("template_id", ""),
                        "tanimoto_to_target": template_exec.get("tanimoto_to_target"),
                    },
                )
            )
        elif not bool(template_exec.get("target_in_products", False)):
            missing_evidence.append(
                _validation_gate_item(
                    "template_target_not_generated",
                    "template_execution",
                    "forward template did not regenerate the target; template coverage may be incomplete",
                    {
                        "template_id": template_exec.get("template_id", ""),
                        "best_match": template_exec.get("best_match", ""),
                        "tanimoto_to_target": template_exec.get("tanimoto_to_target"),
                    },
                )
            )

    scaffold_align = checks.get("scaffold_alignment") or {}
    if scaffold_align and not bool(scaffold_align.get("aligned", True)):
        override_allowed.append(
            _validation_gate_item(
                "scaffold_not_aligned",
                "scaffold_alignment",
                "global MCS scaffold alignment is low; allow override only with site-specific evidence",
                {
                    "coverage_ratio": scaffold_align.get("coverage_ratio"),
                    "summary": scaffold_align.get("summary", ""),
                    "site_audit_pass": (site_audit or {}).get("pass"),
                },
            )
        )

    bond_topo = checks.get("bond_topology") or {}
    bond_topology_violations = [
        item for item in bond_topo.get("violations", []) or []
        if isinstance(item, Mapping)
    ]
    if bond_topo.get("has_hard_fail"):
        hard_blocks.append(
            _validation_gate_item(
                "bond_topology_violation",
                "bond_topology",
                str(bond_topo.get("summary") or "bond topology hard failure"),
                {
                    "violations": [
                        item for item in bond_topology_violations
                        if item.get("severity") == "hard_fail"
                    ] or bond_topology_violations,
                },
            )
        )
    elif bond_topology_violations:
        for item in bond_topology_violations:
            override_allowed.append(
                _validation_gate_item(
                    str(item.get("code") or item.get("type") or "bond_topology_requires_evidence"),
                    "bond_topology",
                    str(item.get("reason") or bond_topo.get("summary") or "bond topology requires atom-source evidence"),
                    {
                        "severity": item.get("severity", ""),
                        "policy_note": bond_topo.get("policy_note", ""),
                    },
                )
            )

    atom_mapping = checks.get("atom_mapping") or {}
    if atom_mapping:
        atom_map_proof_required = [
            item for item in atom_mapping.get("violations", []) or []
            if item.get("severity") == "requires_evidence"
        ]
        atom_map_missing = [
            item for item in atom_mapping.get("violations", []) or []
            if item.get("severity") == "missing_evidence"
        ]
        for item in atom_map_proof_required:
            target_bucket = override_allowed if (
                action_context.get("in_ring")
                or str((ring_topology or {}).get("risk_level", "")) in {"high", "critical"}
            ) else missing_evidence
            target_bucket.append(
                _validation_gate_item(
                    str(item.get("code") or "atom_mapping_requires_evidence"),
                    "atom_mapping",
                    str(item.get("message") or "atom-source mapping requires independent evidence"),
                    item.get("evidence") or {
                        "status": atom_mapping.get("status", ""),
                        "confidence": atom_mapping.get("confidence", ""),
                        "mcs_product_coverage": atom_mapping.get("mcs_product_coverage"),
                    },
                )
            )
        for item in atom_map_missing:
            missing_evidence.append(
                _validation_gate_item(
                    str(item.get("code") or "atom_mapping_missing_evidence"),
                    "atom_mapping",
                    str(item.get("message") or "atom-source mapping evidence is incomplete"),
                    item.get("evidence") or {
                        "status": atom_mapping.get("status", ""),
                        "confidence": atom_mapping.get("confidence", ""),
                    },
                )
            )

    if ring_topology:
        hard_ring_violations = [
            item for item in ring_topology.get("violations", []) or []
            if item.get("severity") == "hard_fail"
        ]
        proof_required = [
            item for item in ring_topology.get("violations", []) or []
            if item.get("severity") == "requires_evidence"
        ]
        for item in hard_ring_violations:
            hard_blocks.append(
                _validation_gate_item(
                    str(item.get("code") or "ring_topology_violation"),
                    "ring_topology",
                    str(item.get("message") or ring_topology.get("summary") or "ring topology hard failure"),
                    {
                        "delta_codes": ring_topology.get("delta_codes", []),
                        "risk_level": ring_topology.get("risk_level", ""),
                    },
                )
            )
        if proof_required:
            template_missed = bool(
                template_exec
                and template_exec.get("attempted")
                and not template_exec.get("target_in_products")
            )
            if action_context.get("in_ring") or template_missed:
                for item in proof_required:
                    override_allowed.append(
                        _validation_gate_item(
                            str(item.get("code") or "ring_topology_requires_evidence"),
                            "ring_topology",
                            str(item.get("message") or "ring topology change requires independent evidence"),
                            {
                                "delta_codes": ring_topology.get("delta_codes", []),
                                "risk_level": ring_topology.get("risk_level", ""),
                                "template_target_in_products": template_exec.get("target_in_products"),
                                "action_in_ring": action_context.get("in_ring"),
                            },
                        )
                    )
            else:
                for item in proof_required:
                    soft_warnings.append(
                        _validation_gate_item(
                            str(item.get("code") or "ring_topology_requires_evidence"),
                            "ring_topology",
                            str(item.get("message") or "ring topology change requires independent evidence"),
                            {
                                "delta_codes": ring_topology.get("delta_codes", []),
                                "risk_level": ring_topology.get("risk_level", ""),
                            },
                        )
                    )

    reaction_specific = checks.get("reaction_specific") or {}
    if reaction_specific.get("has_hard_fail"):
        hard_blocks.append(
            _validation_gate_item(
                "reaction_specific_violation",
                "reaction_specific",
                str(reaction_specific.get("summary") or "reaction-specific hard failure"),
                {"violations": reaction_specific.get("violations", [])},
            )
        )

    ring_family = checks.get("ring_family") or {}
    if ring_family.get("has_hard_fail"):
        for item in ring_family.get("violations", []) or []:
            if item.get("severity") != "hard_fail":
                continue
            hard_blocks.append(
                _validation_gate_item(
                    str(item.get("code") or "ring_family_violation"),
                    "ring_family",
                    str(item.get("message") or ring_family.get("summary") or "ring-family hard failure"),
                    item.get("evidence") or {"violations": ring_family.get("violations", [])},
                )
            )
    elif ring_family.get("violations"):
        for item in ring_family.get("violations", []) or []:
            override_allowed.append(
                _validation_gate_item(
                    str(item.get("code") or "ring_family_requires_evidence"),
                    "ring_family",
                    str(item.get("message") or "ring-family topology evidence is incomplete"),
                    item.get("evidence") or {"violations": ring_family.get("violations", [])},
                )
            )

    fg_compat = checks.get("fg_compatibility") or {}
    fg_warnings = fg_compat.get("warnings", []) or []
    if fg_compat and not bool(fg_compat.get("compatible", True)):
        override_allowed.append(
            _validation_gate_item(
                "forbidden_fg",
                "fg_compatibility",
                "functional-group compatibility is forbidden under the proposed conditions; override requires explicit protection or condition rationale",
                {"warnings": fg_warnings},
            )
        )
    elif fg_warnings:
        soft_warnings.append(
            _validation_gate_item(
                "fg_compatibility_warning",
                "fg_compatibility",
                "functional-group compatibility warnings are present",
                {"warnings": fg_warnings},
            )
        )

    if eas_site.get("applicable"):
        eas_state = str(eas_site.get("state", "") or "")
        eas_evidence = {
            "state": eas_state,
            "reason_code": eas_site.get("reason_code", ""),
            "desired_atom_idx": eas_site.get("desired_atom_idx"),
            "candidate_sites": list(eas_site.get("candidate_sites", []) or []),
            "rule_ids": [
                str(item.get("rule_id", ""))
                for item in eas_site.get("directors", []) or []
                if isinstance(item, Mapping) and item.get("rule_id")
            ],
            "atom_mapping_is_selectivity_proof": False,
        }
        if eas_state == "unsupported":
            override_allowed.append(
                _validation_gate_item(
                    "eas_site_selectivity_requires_evidence",
                    "eas_site_selectivity",
                    str(
                        eas_site.get("summary")
                        or "EAS site selectivity is unsupported or inconclusive."
                    ),
                    eas_evidence,
                    required_evidence=list(
                        eas_site.get("required_evidence", []) or []
                    ),
                )
            )
        elif eas_state == "inconclusive":
            missing_evidence.append(
                _validation_gate_item(
                    "eas_site_selectivity_not_evaluated",
                    "eas_site_selectivity",
                    str(
                        eas_site.get("summary")
                        or "The EAS heuristic does not cover this substrate or site map."
                    ),
                    eas_evidence,
                    required_evidence=list(
                        eas_site.get("required_evidence", []) or []
                    ),
                )
            )
        elif eas_state == "competing":
            soft_warnings.append(
                _validation_gate_item(
                    "eas_site_competition",
                    "eas_site_selectivity",
                    str(
                        eas_site.get("summary")
                        or "Other EAS sites have director support and may compete."
                    ),
                    eas_evidence,
                )
            )

    precursor_state = checks.get("precursor_state") or {}
    if precursor_state.get("status") == "open_shell_detected":
        evidence = {
            "declared": bool(precursor_state.get("declared")),
            "total_radical_electrons": precursor_state.get(
                "total_radical_electrons", 0
            ),
            "open_shell_precursors": precursor_state.get(
                "open_shell_precursors", []
            ),
            "required_evidence": [
                "closed-shell precursor correction or in-situ generation method",
                "lifetime, persistence, or steady-state rationale",
                "atom source and mechanistic role",
                "chemoselectivity and site-selectivity rationale",
            ],
        }
        override_allowed.append(
            _validation_gate_item(
                "open_shell_precursor_requires_evidence",
                "precursor_state",
                "an open-shell carbene, radical, or radical-ion precursor was detected; replace an unintended placeholder with the correct closed-shell precursor, or justify its in-situ generation, lifetime, atom source, mechanistic role, and selectivity",
                evidence,
            )
        )

    policy_adjustments = _apply_controlled_family_downgrades(
        forward_validation,
        hard_blocks=hard_blocks,
        override_allowed=override_allowed,
        missing_evidence=missing_evidence,
    )

    if hard_blocks:
        gate_state = "hard_block"
        commit_policy = "block"
        recommended_action = "Choose another candidate or fix the precursor set; do not override this validation hard block."
    elif override_allowed:
        gate_state = "override_required"
        commit_policy = "requires_override"
        recommended_action = (
            "Treat this as a proof obligation, not automatic disproof: revise the "
            "precursor/action if the chemistry is weak, or commit only with explicit "
            "validation_override grounded in atom-source, tether, anchor, and mechanism evidence."
        )
    elif missing_evidence:
        gate_state = "missing_evidence"
        commit_policy = "allow_with_note"
        recommended_action = "Treat validation as inconclusive; do not use template absence as positive evidence."
    elif soft_warnings:
        gate_state = "warning"
        commit_policy = "allow_with_warning"
        recommended_action = "Commit may proceed only after addressing the warnings in reasoning."
    else:
        gate_state = "pass"
        commit_policy = "allow"
        recommended_action = "Commit may proceed after normal chemistry reasoning and candidate comparison."

    return {
        "gate_state": gate_state,
        "commit_policy": commit_policy,
        "hard_blocks": hard_blocks,
        "soft_warnings": soft_warnings,
        "override_allowed": override_allowed,
        "missing_evidence": missing_evidence,
        "llm_override_allowed": gate_state == "override_required",
        "recommended_action": recommended_action,
        "score": assessment.get("feasibility_score"),
        "legacy_pass": assessment.get("pass"),
        "legacy_hard_fail_reasons": assessment.get("hard_fail_reasons"),
        "policy_adjustments": policy_adjustments,
    }
