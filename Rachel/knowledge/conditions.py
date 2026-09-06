"""Restricted declarative condition evaluation for Rachel knowledge assets."""

from __future__ import annotations

from collections.abc import Collection, Mapping
from typing import Any


class KnowledgeConditionError(ValueError):
    """Raised when a knowledge condition uses unsupported syntax."""


def _field_value(facts: Mapping[str, Any], path: str) -> Any:
    value: Any = facts
    for part in path.split("."):
        if not isinstance(value, Mapping) or part not in value:
            return None
        value = value[part]
    return value


def _values(value: Any) -> list[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Collection):
        return [value]
    return list(value)


def _compare(actual: Any, operation: str, expected: Any) -> bool:
    if operation == "exists":
        return (actual is not None) is bool(expected)
    if operation == "eq":
        return actual == expected
    if operation == "ne":
        return actual != expected
    if operation == "in":
        return actual in _values(expected)
    if operation == "contains":
        return expected in _values(actual)
    if operation == "contains_any":
        return bool(set(_values(actual)) & set(_values(expected)))
    if operation == "contains_all":
        return set(_values(expected)) <= set(_values(actual))
    if operation == "prefix_any":
        prefixes = [str(item) for item in _values(expected)]
        return any(
            str(item).startswith(prefix)
            for item in _values(actual)
            for prefix in prefixes
        )
    if operation in {"lt", "lte", "gt", "gte"}:
        if actual is None or expected is None:
            return False
        if operation == "lt":
            return actual < expected
        if operation == "lte":
            return actual <= expected
        if operation == "gt":
            return actual > expected
        return actual >= expected
    if operation == "range":
        bounds = _values(expected)
        if actual is None or len(bounds) != 2:
            return False
        lower, upper = bounds
        return (lower is None or actual >= lower) and (
            upper is None or actual <= upper
        )
    raise KnowledgeConditionError(
        f"unsupported knowledge condition operation: {operation}"
    )


def condition_matches(condition: Any, facts: Mapping[str, Any]) -> bool:
    """Evaluate a JSON condition without importing or executing pack code."""
    if condition in (None, {}, []):
        return True
    if not isinstance(condition, Mapping):
        raise KnowledgeConditionError("knowledge condition must be an object")

    checks: list[bool] = []
    if "all" in condition:
        clauses = condition["all"]
        if not isinstance(clauses, list):
            raise KnowledgeConditionError("condition all must be a list")
        checks.append(all(condition_matches(item, facts) for item in clauses))
    if "any" in condition:
        clauses = condition["any"]
        if not isinstance(clauses, list):
            raise KnowledgeConditionError("condition any must be a list")
        checks.append(any(condition_matches(item, facts) for item in clauses))
    if "not" in condition:
        checks.append(not condition_matches(condition["not"], facts))

    if "field" in condition:
        field = str(condition.get("field", "")).strip()
        operation = str(condition.get("op", "eq")).strip()
        if not field:
            raise KnowledgeConditionError("condition field cannot be empty")
        checks.append(
            _compare(
                _field_value(facts, field),
                operation,
                condition.get("value"),
            )
        )

    shorthand: list[bool] = []
    tags = set(_values(facts.get("tags", [])))
    events = {str(item) for item in _values(facts.get("events", []))}
    if "any_tags" in condition:
        shorthand.append(bool(tags & set(_values(condition["any_tags"]))))
    if "all_tags" in condition:
        shorthand.append(set(_values(condition["all_tags"])) <= tags)
    if "any_events" in condition:
        shorthand.append(bool(events & set(_values(condition["any_events"]))))
    if "all_events" in condition:
        shorthand.append(set(_values(condition["all_events"])) <= events)
    if "event_prefixes" in condition:
        prefixes = [str(item) for item in _values(condition["event_prefixes"])]
        shorthand.append(
            any(event.startswith(prefix) for event in events for prefix in prefixes)
        )
    if shorthand:
        checks.append(any(shorthand))

    if not checks:
        unknown = sorted(set(condition) - {"description"})
        if unknown:
            raise KnowledgeConditionError(
                f"unsupported knowledge condition keys: {', '.join(unknown)}"
            )
        return True
    return all(checks)


__all__ = [
    "KnowledgeConditionError",
    "condition_matches",
]
