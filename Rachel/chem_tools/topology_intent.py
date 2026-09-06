"""Structured topology-intent detection shared by validation components."""

from __future__ import annotations

from typing import Any, Mapping, Set


_TOPOLOGY_DELTAS = {
    "ring_closure",
    "ring_opening",
    "ring_construction",
    "ring_system_merge",
    "new_fused_ring_system",
    "new_spiro_center",
    "new_bridged_system",
    "new_medium_ring",
    "new_macrocycle",
    "scaffold_rewrite",
}

_TOPOLOGY_RISK_TAGS = {
    "ring_construction",
    "ring_opening",
    "scaffold_edit",
    "topology_change",
}

_NO_CHANGE_VALUES = {"none", "no_change", "preserve", "preserved"}


def _normalized_values(raw: Any) -> Set[str]:
    if not isinstance(raw, (list, tuple, set)):
        return set()
    return {
        str(item or "").strip().lower()
        for item in raw
        if str(item or "").strip()
    }


def action_declares_topology_change(action_context: Mapping[str, Any]) -> bool:
    """Use explicit deltas/risk fields, never reaction-name substrings."""
    expected = str(action_context.get("expected_ring_change", "") or "").strip().lower()
    if expected and expected not in _NO_CHANGE_VALUES:
        return True
    if _normalized_values(action_context.get("intended_deltas")) & _TOPOLOGY_DELTAS:
        return True
    return bool(_normalized_values(action_context.get("risk_tags")) & _TOPOLOGY_RISK_TAGS)
