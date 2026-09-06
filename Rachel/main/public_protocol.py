from __future__ import annotations

from typing import Any, Dict


PUBLIC_FIELD_KEY_MAP = {
    "candidate_id": "action_id",
    "candidate_ids": "action_ids",
    "candidate_summary": "action_summary",
    "candidate_comparison": "action_comparison",
    "selected_candidate_id": "selected_action_id",
    "rejected_candidates": "rejected_actions",
    "rejected_candidate_ids": "rejected_action_ids",
    "rescue_steps": "continuation_steps",
    "rescue_id": "continuation_id",
    "rescue_step_idx": "continuation_step_idx",
    "rescue_step_count": "continuation_step_count",
    "rescue_continuation": "strategy_continuation",
    "rescue_continuation_brief": "strategy_continuation_brief",
    "rescue_controls": "continuation_controls",
    "rescue_status": "continuation_status",
    "rescue_abort": "continuation_abort",
    "active_rescue_continuations": "active_strategy_continuations",
    "blocked_rescue_count": "pending_continuation_count",
    "total_rescue_continuations": "total_strategy_continuations",
}

PUBLIC_LITERAL_MAP = {
    "strategy.rescue_continuation_active": "strategy.continuation_active",
    "route_sketch_multi_step_rescue": "route_sketch_multi_step_continuation",
    "rescue_continuation": "strategy_continuation",
    "rescue_continuation_pending": "strategy_continuation_pending",
    "rescue_continuation_not_found": "strategy_continuation_not_found",
    "rescue_continuation_unbound": "strategy_continuation_unbound",
    "terminal_rescue_continuation_required": "terminal_strategy_continuation_required",
    "rescue_aborted": "strategy_continuation_aborted",
    "rescue-abort": "strategy-continuation-abort",
    "commit credible attempt to create rescue_continuation, or force_accept_without_rescue with reason": "commit credible attempt to create strategy_continuation, or force_accept_without_rescue with reason",
}


def _project_public_string(value: str) -> str:
    if value.startswith("rescue:"):
        return f"continuation:{value[len('rescue:') :]}"
    if value in PUBLIC_LITERAL_MAP:
        return PUBLIC_LITERAL_MAP[value]
    replacements = (
        ("rescue_abort(rescue_id, reason)", "continuation_abort(continuation_id, reason)"),
        ("rescue_abort(reason=...)", "continuation_abort(reason=...)"),
        ("rescue_status()", "continuation_status()"),
    )
    projected = value
    for old, new in replacements:
        projected = projected.replace(old, new)
    return projected


def internal_continuation_id(value: Any) -> str:
    """Map a public continuation id back to the legacy internal id."""
    text = str(value or "").strip()
    if text.startswith("continuation:"):
        return f"rescue:{text[len('continuation:') :]}"
    return text


def project_public_payload(value: Any) -> Any:
    """Project legacy internal fields into the current public protocol."""
    if isinstance(value, str):
        return _project_public_string(value)
    if isinstance(value, list):
        return [project_public_payload(item) for item in value]
    if not isinstance(value, dict):
        return value

    projected: Dict[str, Any] = {}
    for key, raw_item in value.items():
        public_key = PUBLIC_FIELD_KEY_MAP.get(key, key)
        item = project_public_payload(raw_item)
        if public_key in projected:
            if key == public_key:
                projected[public_key] = item
            elif projected[public_key] in ("", [], {}, None):
                projected[public_key] = item
            continue
        projected[public_key] = item
    return projected
