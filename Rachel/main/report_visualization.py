"""Lightweight human-report annotations for visualization.json."""

from __future__ import annotations

from typing import Any, Dict


def annotate_visualization(
    payload: Dict[str, Any],
    report_view: Dict[str, Any],
) -> Dict[str, Any]:
    steps = {
        str(step.get("step_id") or ""): step
        for step in report_view.get("steps", [])
    }
    for node in payload.get("nodes", []):
        if node.get("type") != "reaction":
            continue
        step = steps.get(str(node.get("canonical_id") or node.get("id") or ""))
        if not step:
            continue
        conditions = step.get("conditions") or {}
        node["forward_step"] = step.get("forward_step")
        node["validation_gate_state"] = (
            (step.get("validation") or {}).get("state") or "not_recorded"
        )
        node["process_condition_coverage"] = {
            "specified": int(conditions.get("specified_count", 0) or 0),
            "required": int(conditions.get("required_count", 0) or 0),
            "planning_suggestions": len(
                conditions.get("planning_suggestions", []) or []
            ),
            "unknown": len(conditions.get("unknown_fields", []) or []),
            "missing": len(conditions.get("missing_fields", []) or []),
            "invalid": len(conditions.get("invalid_fields", []) or []),
        }
    payload.setdefault("meta", {})["forward_step_count"] = int(
        report_view.get("step_count", 0) or 0
    )
    return payload
