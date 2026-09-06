"""Shared validation finding schema.

Audit modules should report facts as findings. Policy code can later decide
whether those findings become hard blocks, overrides, warnings, or missing
evidence in a specific reaction context.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Optional


VALID_SEVERITIES = frozenset({
    "info",
    "warning",
    "missing_evidence",
    "requires_evidence",
    "hard_fail",
})


def make_finding(
    *,
    code: str,
    severity: str,
    source: str,
    message: str,
    evidence: Optional[Mapping[str, Any]] = None,
    required_evidence: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Create a normalized validation finding.

    The returned dict intentionally uses plain JSON-compatible containers
    because findings are stored in session artifacts and reports.
    """
    normalized_severity = severity if severity in VALID_SEVERITIES else "warning"
    return {
        "code": str(code or "validation_finding"),
        "severity": normalized_severity,
        "source": str(source or "validation"),
        "message": str(message or ""),
        "evidence": dict(evidence or {}),
        "required_evidence": list(required_evidence or []),
    }

