"""Deterministic contract for LLM-declared route-step risk information."""

from __future__ import annotations

import json
from typing import Any, Dict, List


PROCESS_CONDITION_FIELDS = (
    "solvent",
    "equivalents",
    "addition_order",
    "temperature",
    "atmosphere",
    "reaction_time",
    "workup",
)

OPTIONAL_PROCESS_CONDITION_FIELDS = (
    "scale",
    "concentration",
    "pressure",
)

_RECOGNIZED_PROCESS_CONDITION_FIELDS = (
    "scale",
    "concentration",
    "solvent",
    "equivalents",
    "addition_order",
    "temperature",
    "pressure",
    "atmosphere",
    "reaction_time",
    "workup",
)

PROCESS_CONDITIONS_INSUFFICIENT = (
    "conditions_insufficient_for_process_assessment"
)
PROCESS_CONDITIONS_PRELIMINARY_REVIEW = (
    "conditions_available_for_preliminary_process_review"
)

PROCESS_CONDITION_BASIS_KINDS = (
    "user_provided",
    "literature",
    "experimental",
    "sop",
    "planning_hypothesis",
)
EVIDENCE_BACKED_BASIS_KINDS = {
    "user_provided",
    "literature",
    "experimental",
    "sop",
}

_NON_VALUES = {
    "",
    "n/a",
    "na",
    "not known",
    "not provided",
    "not specified",
    "tbd",
    "to be determined",
    "unknown",
    "unspecified",
}


def _has_concrete_value(value: Any) -> bool:
    if value is None or isinstance(value, bool):
        return False
    if isinstance(value, str):
        return value.strip().lower() not in _NON_VALUES
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, dict):
        return bool(value) and all(
            str(key).strip() and _has_concrete_value(item)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple)):
        return bool(value) and all(_has_concrete_value(item) for item in value)
    return False


def assess_process_conditions(conditions: Any) -> Dict[str, Any]:
    """Check structured field coverage without judging chemical correctness."""
    declared = conditions if isinstance(conditions, dict) else {}
    specified: List[str] = []
    core_specified: List[str] = []
    optional_specified: List[str] = []
    evidence_backed: List[str] = []
    planning_suggestions: List[str] = []
    unclassified_specified: List[str] = []
    unknown: List[str] = []
    optional_unknown: List[str] = []
    missing: List[str] = []
    invalid: List[str] = []
    optional_invalid: List[str] = []

    for field_name in _RECOGNIZED_PROCESS_CONDITION_FIELDS:
        is_required = field_name in PROCESS_CONDITION_FIELDS
        if field_name not in declared:
            if is_required:
                missing.append(field_name)
            continue
        entry = declared.get(field_name)
        if not isinstance(entry, dict):
            (invalid if is_required else optional_invalid).append(field_name)
            continue
        status = str(entry.get("status") or "").strip().lower()
        if status == "unknown" and _has_concrete_value(entry.get("reason")):
            (unknown if is_required else optional_unknown).append(field_name)
        elif status == "specified" and _has_concrete_value(entry.get("value")):
            specified.append(field_name)
            (core_specified if is_required else optional_specified).append(field_name)
            basis_kind = str(entry.get("basis_kind") or "").strip().lower()
            if basis_kind in EVIDENCE_BACKED_BASIS_KINDS:
                evidence_backed.append(field_name)
            elif basis_kind == "planning_hypothesis":
                planning_suggestions.append(field_name)
            else:
                # Existing sessions predate basis_kind. Preserve their coverage
                # while making the unclassified provenance visible to reports.
                unclassified_specified.append(field_name)
        else:
            (invalid if is_required else optional_invalid).append(field_name)

    ready = len(core_specified) == len(PROCESS_CONDITION_FIELDS)
    return {
        "required_fields": list(PROCESS_CONDITION_FIELDS),
        "optional_fields": list(OPTIONAL_PROCESS_CONDITION_FIELDS),
        "specified_fields": specified,
        "core_specified_fields": core_specified,
        "optional_specified_fields": optional_specified,
        "evidence_backed_fields": evidence_backed,
        "planning_suggestion_fields": planning_suggestions,
        "unclassified_specified_fields": unclassified_specified,
        "unknown_fields": unknown,
        "optional_unknown_fields": optional_unknown,
        "missing_fields": missing,
        "invalid_fields": invalid,
        "optional_invalid_fields": optional_invalid,
        "additional_fields": sorted(
            str(key)
            for key in declared
            if key not in _RECOGNIZED_PROCESS_CONDITION_FIELDS
        ),
        "specified_count": len(core_specified),
        "declared_specified_count": len(specified),
        "optional_specified_count": len(optional_specified),
        "evidence_backed_count": len(evidence_backed),
        "planning_suggestion_count": len(planning_suggestions),
        "unclassified_specified_count": len(unclassified_specified),
        "required_count": len(PROCESS_CONDITION_FIELDS),
        "optional_count": len(OPTIONAL_PROCESS_CONDITION_FIELDS),
        "coverage_fraction": round(
            len(core_specified) / len(PROCESS_CONDITION_FIELDS), 3
        ),
        "process_assessment_status": (
            PROCESS_CONDITIONS_PRELIMINARY_REVIEW
            if ready
            else PROCESS_CONDITIONS_INSUFFICIENT
        ),
        "semantic_validation": "not_performed",
    }


def format_process_condition_gaps(coverage: Dict[str, Any]) -> str:
    parts = []
    for field_name, label in (
        ("missing_fields", "missing"),
        ("unknown_fields", "unknown"),
        ("invalid_fields", "invalid"),
    ):
        values = coverage.get(field_name) or []
        if values:
            parts.append(f"{label}={','.join(str(item) for item in values)}")
    return "; ".join(parts) or "none"


def format_process_conditions(conditions: Any) -> str:
    declared = conditions if isinstance(conditions, dict) else {}
    rendered = []
    for field_name in _RECOGNIZED_PROCESS_CONDITION_FIELDS:
        if (
            field_name in OPTIONAL_PROCESS_CONDITION_FIELDS
            and field_name not in declared
        ):
            continue
        entry = declared.get(field_name)
        if not isinstance(entry, dict):
            rendered.append(f"{field_name}=missing")
            continue
        status = str(entry.get("status") or "").strip().lower()
        if status == "specified" and _has_concrete_value(entry.get("value")):
            value = json.dumps(entry.get("value"), ensure_ascii=False, sort_keys=True)
            rendered.append(f"{field_name}={value}")
        elif status == "unknown":
            rendered.append(f"{field_name}=unknown")
        else:
            rendered.append(f"{field_name}=invalid")
    return "; ".join(rendered)
