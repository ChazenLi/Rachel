"""Human-oriented route report projection and renderers.

The session, route tree, and post-route audit remain the complete audit record.
This module builds one smaller forward-route view for chemists and renders that
same view as text, Markdown, and HTML.
"""

from __future__ import annotations

import html
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from Rachel.chem_tools.validation_contract import build_validation_contract

from .retro_tree import ReactionNode, RetrosynthesisTree
from .risk_protocol import assess_process_conditions


def _node_id_sort_key(node_id: str) -> tuple[int, str]:
    try:
        return (int(str(node_id).rsplit("_", 1)[-1]), str(node_id))
    except (TypeError, ValueError):
        return (10**9, str(node_id))


def forward_reactions(tree: RetrosynthesisTree) -> List[ReactionNode]:
    product_to_rxn = {rxn.product_node: rxn for rxn in tree.reaction_nodes}
    by_id = {rxn.step_id: rxn for rxn in tree.reaction_nodes}
    in_degree = {rxn.step_id: 0 for rxn in tree.reaction_nodes}
    forward: Dict[str, List[str]] = defaultdict(list)
    for rxn in tree.reaction_nodes:
        for reactant_id in rxn.reactant_nodes:
            dependency = product_to_rxn.get(reactant_id)
            if dependency is not None:
                in_degree[rxn.step_id] += 1
                forward[dependency.step_id].append(rxn.step_id)

    queue = sorted(step_id for step_id, degree in in_degree.items() if degree == 0)
    ordered: List[ReactionNode] = []
    while queue:
        step_id = queue.pop(0)
        ordered.append(by_id[step_id])
        for neighbor in sorted(forward[step_id]):
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)
                queue.sort()
    return ordered


def _as_dict(item: Any) -> Dict[str, Any]:
    if isinstance(item, dict):
        return dict(item)
    to_dict = getattr(item, "to_dict", None)
    return dict(to_dict()) if callable(to_dict) else {}


def _format_rejected(item: Any) -> str:
    if isinstance(item, str):
        return item
    if isinstance(item, dict):
        name = (
            item.get("method")
            or item.get("action_id")
            or item.get("candidate_id")
            or item.get("reaction_type")
            or item.get("name")
            or "alternative"
        )
        reason = item.get("reason") or item.get("rationale") or item.get("note")
        return f"{name}: {reason}" if reason else str(name)
    return str(item)


def _route_plan_view(
    tree: RetrosynthesisTree,
    route_plan: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    if isinstance(route_plan, dict) and route_plan.get("route_thesis"):
        return dict(route_plan)

    selected: Dict[str, Any] = {}
    selected_revision = -1
    for rxn in tree.reaction_nodes:
        decision = rxn.llm_decision
        audit = decision.decision_audit if decision else {}
        candidate = (
            (audit or {}).get("route_plan_snapshot")
            or (audit or {}).get("route_plan_brief")
            or {}
        )
        if not isinstance(candidate, dict) or not candidate.get("route_thesis"):
            continue
        try:
            revision = int(candidate.get("revision", 0) or 0)
        except (TypeError, ValueError):
            revision = 0
        if revision >= selected_revision:
            selected = dict(candidate)
            selected_revision = revision
    return selected


def _finding_codes(items: Iterable[Any]) -> List[str]:
    codes: List[str] = []
    for item in items or []:
        if isinstance(item, dict):
            code = item.get("code") or item.get("reason") or item.get("name")
        else:
            code = item
        text = str(code or "").strip()
        if text and text not in codes:
            codes.append(text)
    return codes


def _validation_view(rxn: ReactionNode) -> Dict[str, Any]:
    decision = rxn.llm_decision
    audit = decision.decision_audit if decision else {}
    gate = (audit or {}).get("validation_gate") or {}
    if not gate and isinstance(rxn.forward_validation, dict):
        gate = rxn.forward_validation.get("gate") or {}
    if not gate:
        return {"state": "not_recorded", "key_findings": []}

    contract = build_validation_contract({"gate": gate})
    findings: List[Dict[str, str]] = []
    for category in (
        "contradictions",
        "proof_obligations",
        "evidence_gaps",
        "tool_limits",
        "warnings",
        "system_errors",
    ):
        for code in _finding_codes(contract.get(category, [])):
            findings.append({"category": category, "code": code})
    return {
        "state": (contract.get("decision_gate") or {}).get("state", "not_recorded"),
        "key_findings": findings[:6],
    }


def _condition_values(
    conditions: Dict[str, Any],
    fields: Iterable[str],
) -> List[Dict[str, Any]]:
    values = []
    for field_name in fields:
        entry = conditions.get(field_name) or {}
        values.append({
            "field": field_name,
            "value": entry.get("value"),
            "basis_kind": str(entry.get("basis_kind") or "").strip().lower(),
            "basis": str(entry.get("basis") or "").strip(),
        })
    return values


def _condition_view(conditions: Any) -> Dict[str, Any]:
    declared = conditions if isinstance(conditions, dict) else {}
    coverage = assess_process_conditions(declared)
    reasons: List[str] = []
    for field_name in coverage.get("unknown_fields", []):
        reason = str((declared.get(field_name) or {}).get("reason") or "").strip()
        if reason and reason not in reasons:
            reasons.append(reason)
    return {
        "status": coverage["process_assessment_status"],
        "specified_count": coverage["specified_count"],
        "required_count": coverage["required_count"],
        "evidence_backed": _condition_values(
            declared, coverage.get("evidence_backed_fields", [])
        ),
        "planning_suggestions": _condition_values(
            declared, coverage.get("planning_suggestion_fields", [])
        ),
        "unclassified_specified": _condition_values(
            declared, coverage.get("unclassified_specified_fields", [])
        ),
        "unknown_fields": list(coverage.get("unknown_fields", [])),
        "unknown_reasons": reasons,
        "missing_fields": list(coverage.get("missing_fields", [])),
        "invalid_fields": list(coverage.get("invalid_fields", [])),
    }


def _frontier_view(payload: Dict[str, Any]) -> Optional[Dict[str, str]]:
    if payload.get("hypothesis_status") != "unprecedented_hypothesis":
        return None
    return {
        "status": "unprecedented_hypothesis",
        "basis": str(payload.get("hypothesis_basis") or "").strip(),
    }


def _scouting_rows(scouting_view: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    reviews = (scouting_view or {}).get("latest_review_by_node", {}) or {}
    counts = (scouting_view or {}).get("round_count_by_node", {}) or {}
    seeds = (scouting_view or {}).get("deferred_review_seeds", {}) or {}
    adopted_counts: Counter[str] = Counter()
    for item in (scouting_view or {}).get("adopted_sources", []) or []:
        adopted_counts[str(item.get("node_id") or "")] += 1

    rows = []
    for node_id, review in reviews.items():
        shortlist = []
        for item in review.get("shortlisted_directions", []) or []:
            direction = item.get("primary_direction") or {}
            shortlist.append({
                "summary": str(direction.get("direction_summary") or "").strip(),
                "frontier": _frontier_view(direction),
            })
        deferred = []
        for item in seeds.get(node_id, []) or []:
            deferred.append({
                "summary": str(item.get("direction_summary") or "").strip(),
                "frontier": _frontier_view(item),
            })
        rows.append({
            "node_id": str(node_id),
            "round_count": int(counts.get(node_id, 0) or 0),
            "review_summary": str(review.get("review_summary") or "").strip(),
            "shortlist": shortlist,
            "deferred_review_seeds": deferred,
            "adopted_count": adopted_counts[str(node_id)],
        })
    return rows


def _audit_summary(payload: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not payload:
        return {}
    buyability = payload.get("buyability") or {}
    safety = payload.get("safety") or {}
    safety_summary = safety.get("summary") or {}
    grades = (buyability.get("summary") or {}).get("by_evidence_grade") or {}
    return {
        "status": str(payload.get("status") or "unknown"),
        "incompleteness_reasons": list(payload.get("incompleteness_reasons") or []),
        "buyability_grades": {
            grade: int(grades.get(grade, 0) or 0) for grade in ("B0", "B1", "B2")
        },
        "s1_substances": int(
            safety_summary.get("substances_with_source_attributed_ghs", 0) or 0
        ),
        "steps_missing_risk": int(
            safety_summary.get("steps_missing_explicit_risk_assessment", 0) or 0
        ),
        "conditions_ready_steps": int(
            safety_summary.get(
                "steps_with_process_conditions_for_preliminary_review", 0
            ) or 0
        ),
        "conditions_incomplete_steps": int(
            safety_summary.get("steps_with_incomplete_process_conditions", 0) or 0
        ),
        "markdown_file": "post_route_audit.md",
        "json_file": "post_route_audit.json",
    }


def _decision_reason(record: Dict[str, Any]) -> str:
    detail = str(record.get("reasoning_detail") or "").strip()
    if detail:
        return detail
    summary = str(record.get("reasoning_summary") or "").strip()
    if len(summary) >= 120:
        summary = summary.rsplit(" ", 1)[0] or summary
        return summary.rstrip(".,;:") + "..."
    return summary


def _decision_report_views(
    decision_history: Optional[Iterable[Any]],
) -> tuple[Dict[str, Dict[str, str]], List[Dict[str, str]]]:
    terminal_bases: Dict[str, Dict[str, str]] = {}
    unresolved_skips: Dict[str, Dict[str, str]] = {}
    terminal_reviews: List[Dict[str, str]] = []
    terminal_labels = {
        "accept-terminal": "Accepted as a planned starting material",
        "skip": "Planning stopped at this molecule",
        "rescue-abort": "Terminal after continuation review",
    }
    for item in decision_history or []:
        record = _as_dict(item)
        action = str(record.get("action") or "").strip()
        molecule = str(record.get("molecule") or "").strip()
        if not molecule:
            continue
        if action == "skip":
            basis = {
                "source": terminal_labels[action],
                "reason": _decision_reason(record),
                "status": str(record.get("outcome") or "skipped"),
                "action": action,
                "kind": "planning_stop",
            }
            unresolved_skips[molecule] = basis
            terminal_bases[molecule] = basis
            continue
        if action in {"accept-terminal", "commit"}:
            unresolved_skips.pop(molecule, None)
            if action == "commit":
                terminal_bases.pop(molecule, None)
                continue
        if action == "review-terminal":
            terminal_bases.pop(molecule, None)
            terminal_reviews.append({
                "kind": "terminal",
                "title": "Terminal reopened for further planning",
                "smiles": molecule,
                "reason": _decision_reason(record),
                "status": "reopened",
            })
            continue
        if action == "rescue-abort" and molecule in unresolved_skips:
            terminal_bases[molecule] = unresolved_skips[molecule]
            continue
        if action in terminal_labels:
            terminal_bases[molecule] = {
                "source": terminal_labels[action],
                "reason": _decision_reason(record),
                "status": str(record.get("outcome") or "terminal"),
                "action": action,
                "kind": "terminal",
            }
    return terminal_bases, terminal_reviews


def _variant_review_view(variant_info: Optional[Dict[str, Any]]) -> Dict[str, str]:
    if not isinstance(variant_info, dict) or not variant_info:
        return {}
    mode = str(variant_info.get("review_mode") or "").strip()
    mode_labels = {
        "replace_subtree": "Replace reviewed subtree",
        "extend_terminal": "Extend reviewed terminal",
    }
    return {
        "kind": "node",
        "title": "Node review variant",
        "node_id": str(variant_info.get("reviewed_node_id") or "").strip(),
        "smiles": str(variant_info.get("reviewed_smiles") or "").strip(),
        "review_mode": mode,
        "review_mode_label": mode_labels.get(mode, mode.replace("_", " ")),
        "reason": str(variant_info.get("reason") or "").strip(),
        "instruction": str(variant_info.get("instruction_summary") or "").strip(),
    }


def build_route_report_view(
    tree: RetrosynthesisTree,
    *,
    route_plan: Optional[Dict[str, Any]] = None,
    scouting_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
    experimental_outcomes: Optional[Iterable[Any]] = None,
    variant_info: Optional[Dict[str, Any]] = None,
    decision_history: Optional[Iterable[Any]] = None,
) -> Dict[str, Any]:
    """Build the single human-report projection from complete route state."""
    outcomes = [_as_dict(item) for item in (experimental_outcomes or [])]
    outcomes = [item for item in outcomes if item]
    outcomes_by_step: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for outcome in outcomes:
        outcomes_by_step[str(outcome.get("step_id") or "")].append(outcome)

    terminal_bases, terminal_reviews = _decision_report_views(decision_history)
    terminals = []
    planning_stops = []
    for node in tree.get_starting_material_nodes():
        terminal_basis = dict(terminal_bases.get(node.smiles) or {})
        item = {
            "node_id": node.node_id,
            "smiles": node.smiles,
            "cs_score": node.cs_score,
            "classification": (
                (node.complexity or {}).get("classification", "")
            ),
            "terminal_basis": terminal_basis,
        }
        if terminal_basis.get("kind") == "planning_stop":
            planning_stops.append(item)
        else:
            terminals.append(item)
    terminals.sort(key=lambda item: _node_id_sort_key(item["node_id"]))
    planning_stops.sort(key=lambda item: _node_id_sort_key(item["node_id"]))

    steps = []
    for forward_step, rxn in enumerate(forward_reactions(tree), 1):
        decision = rxn.llm_decision
        audit = decision.decision_audit if decision else {}
        reactants = []
        for node_id in rxn.reactant_nodes:
            node = tree.molecule_nodes.get(node_id)
            reactants.append({
                "node_id": node_id,
                "smiles": node.smiles if node else node_id,
                "role": node.role if node else "unknown",
            })
        product = tree.molecule_nodes.get(rxn.product_node)
        rejected = decision.rejected_alternatives if decision else []
        if isinstance(rejected, str):
            rejected = [rejected]
        steps.append({
            "forward_step": forward_step,
            "step_id": rxn.step_id,
            "reaction_name": (
                rxn.reaction_type
                or (rxn.template_evidence.template_name if rxn.template_evidence else "")
                or "Unspecified reaction"
            ),
            "reaction_smiles": rxn.reaction_smiles,
            "reactants": reactants,
            "reagents": list(rxn.reagents),
            "product": {
                "node_id": rxn.product_node,
                "smiles": product.smiles if product else rxn.product_node,
            },
            "reasoning": decision.selection_reasoning if decision else "",
            "risk_assessment": (
                decision.risk_assessment
                if decision and decision.risk_assessment
                else (
                    "No explicit assessment was recorded; verify exact conditions, "
                    "SDS, SOP, and scale-up risks before execution"
                )
            ),
            "conditions": _condition_view(
                getattr(decision, "process_conditions", {}) if decision else {}
            ),
            "validation": _validation_view(rxn),
            "validation_override": dict((audit or {}).get("validation_override") or {}),
            "rejected_alternatives": [_format_rejected(item) for item in rejected],
            "experimental_outcomes": outcomes_by_step.get(rxn.step_id, []),
        })

    gate_counts = Counter(step["validation"]["state"] for step in steps)
    outcome_counts = Counter(str(item.get("outcome") or "unknown") for item in outcomes)
    variant_review = _variant_review_view(variant_info)
    route_reviews = ([variant_review] if variant_review else []) + terminal_reviews
    return {
        "schema": "rachel.route_report_view.v1",
        "target": tree.target,
        "target_name": tree.target_name or tree.target,
        "route_status": str(tree.status),
        "route_status_statement": (
            "Route planning is incomplete because unresolved planning stops remain."
            if planning_stops
            else (
                "Route construction complete."
                if str(tree.status) == "complete"
                else "Route planning in progress."
            )
        ),
        "step_count": len(steps),
        "starting_material_count": len(terminals),
        "planning_stop_count": len(planning_stops),
        "route_plan": _route_plan_view(tree, route_plan),
        "route_summary": str(tree.llm_summary or "").strip(),
        "validation_gate_counts": dict(sorted(gate_counts.items())),
        "experimental_outcome_counts": dict(sorted(outcome_counts.items())),
        "steps_without_experimental_outcome": sum(
            1 for step in steps if not step["experimental_outcomes"]
        ),
        "starting_materials": terminals,
        "planning_stops": planning_stops,
        "steps": steps,
        "route_reviews": route_reviews,
        "scouting_reviews": _scouting_rows(scouting_view),
        "post_route_audit": _audit_summary(post_route_audit),
    }


_CONDITION_FIELD_LABELS = {
    "scale": "Scale",
    "concentration": "Concentration",
    "solvent": "Solvent",
    "equivalents": "Equivalents",
    "addition_order": "Addition and operation sequence",
    "temperature": "Temperature",
    "pressure": "Pressure",
    "atmosphere": "Atmosphere",
    "reaction_time": "Reaction time",
    "workup": "Workup",
}

_WORKUP_LABELS = {
    "quench": "Quench",
    "isolation": "Isolation",
    "purification": "Purification",
}

_CONDITION_GROUPS = (
    ("Evidence-backed", "evidence_backed"),
    ("Planning suggestions (unverified)", "planning_suggestions"),
    ("Source unclassified", "unclassified_specified"),
)


def _condition_label(field_name: Any) -> str:
    text = str(field_name or "").strip()
    return _CONDITION_FIELD_LABELS.get(text, text.replace("_", " "))


def _condition_value_text(value: Any) -> str:
    if value is None:
        return "Not specified"
    if isinstance(value, dict):
        return "; ".join(
            f"{_condition_label(key)}: {_condition_value_text(item)}"
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple)):
        return "; ".join(_condition_value_text(item) for item in value)
    return str(value)


def _condition_item_view(item: Dict[str, Any]) -> Dict[str, Any]:
    field_name = str(item.get("field") or "")
    value = item.get("value")
    rendered = {
        "label": _condition_label(field_name),
        "kind": "value",
        "value": _condition_value_text(value),
    }
    if field_name == "addition_order" and isinstance(value, (list, tuple)):
        rendered.update({
            "kind": "steps",
            "items": [_condition_value_text(entry) for entry in value],
        })
    elif field_name == "equivalents" and isinstance(value, dict):
        rendered.update({
            "kind": "pairs",
            "items": [
                {
                    "label": str(key).replace("_", " "),
                    "value": _condition_value_text(entry),
                }
                for key, entry in value.items()
            ],
        })
    elif field_name == "workup" and isinstance(value, dict):
        ordered_keys = [key for key in _WORKUP_LABELS if key in value]
        ordered_keys.extend(key for key in value if key not in _WORKUP_LABELS)
        rendered.update({
            "kind": "pairs",
            "items": [
                {
                    "label": _WORKUP_LABELS.get(key, _condition_label(key)),
                    "value": _condition_value_text(value[key]),
                }
                for key in ordered_keys
            ],
        })
    return rendered


def _condition_display_view(conditions: Dict[str, Any]) -> Dict[str, Any]:
    required = int(conditions.get("required_count", 0) or 0)
    specified = int(conditions.get("specified_count", 0) or 0)
    groups = []
    for label, key in _CONDITION_GROUPS:
        items = conditions.get(key, []) or []
        if items:
            groups.append({
                "label": label,
                "items": [_condition_item_view(item) for item in items],
            })

    issues = []
    for key, label in (
        ("unknown_fields", "Unknown fields"),
        ("missing_fields", "Missing fields"),
        ("invalid_fields", "Invalid fields"),
    ):
        fields = conditions.get(key, []) or []
        if fields:
            issues.append({
                "label": label,
                "value": (
                    f"{len(fields)}/{required} ("
                    + ", ".join(_condition_label(field) for field in fields)
                    + ")"
                ),
            })
    reasons = [str(reason) for reason in (conditions.get("unknown_reasons", []) or [])]
    return {
        "coverage": f"{specified}/{required} specified",
        "groups": groups,
        "issues": issues,
        "unknown_reasons": reasons[:2],
    }


def _append_text_condition_item(lines: List[str], item: Dict[str, Any]) -> None:
    if item["kind"] == "value":
        lines.append(f"      {item['label']}: {item['value']}")
        return
    lines.append(f"      {item['label']}:")
    marker = "{index}." if item["kind"] == "steps" else "-"
    for index, entry in enumerate(item["items"], 1):
        if isinstance(entry, dict):
            value = f"{entry['label']}: {entry['value']}"
        else:
            value = str(entry)
        lines.append(f"        {marker.format(index=index)} {value}")


def _text_condition_lines(conditions: Dict[str, Any]) -> List[str]:
    display = _condition_display_view(conditions)
    lines = [f"    Core fields: {display['coverage']}"]
    for group in display["groups"]:
        lines.append(f"    {group['label']}:")
        for item in group["items"]:
            _append_text_condition_item(lines, item)
    for issue in display["issues"]:
        lines.append(f"    {issue['label']}: {issue['value']}")
    for reason in display["unknown_reasons"]:
        lines.append(f"    Unknown reason: {reason}")
    return lines


def _markdown_condition_lines(conditions: Dict[str, Any]) -> List[str]:
    display = _condition_display_view(conditions)
    lines = ["- **Process conditions**:", f"  - **Core fields**: {display['coverage']}"]
    for group in display["groups"]:
        lines.append(f"  - **{group['label']}**:")
        for item in group["items"]:
            if item["kind"] == "value":
                lines.append(f"    - **{item['label']}**: {item['value']}")
                continue
            lines.append(f"    - **{item['label']}**:")
            for index, entry in enumerate(item["items"], 1):
                if isinstance(entry, dict):
                    lines.append(f"      - **{entry['label']}**: {entry['value']}")
                else:
                    lines.append(f"      {index}. {entry}")
    for issue in display["issues"]:
        lines.append(f"  - **{issue['label']}**: {issue['value']}")
    for reason in display["unknown_reasons"]:
        lines.append(f"  - **Unknown reason**: {reason}")
    return lines


def _html_condition_block(conditions: Dict[str, Any]) -> str:
    display = _condition_display_view(conditions)
    chunks = [
        "<div class='conditions detail-panel'>",
        "<p><strong>Core fields</strong>: ",
        html.escape(display["coverage"]),
        "</p>",
    ]
    for group in display["groups"]:
        chunks.append(f"<h4>{html.escape(group['label'])}</h4><dl>")
        for item in group["items"]:
            chunks.append(f"<dt>{html.escape(item['label'])}</dt><dd>")
            if item["kind"] == "value":
                chunks.append(html.escape(item["value"]))
            elif item["kind"] == "steps":
                chunks.append("<ol>")
                chunks.extend(f"<li>{html.escape(entry)}</li>" for entry in item["items"])
                chunks.append("</ol>")
            else:
                chunks.append("<ul>")
                chunks.extend(
                    "<li><strong>"
                    + html.escape(entry["label"])
                    + "</strong>: "
                    + html.escape(entry["value"])
                    + "</li>"
                    for entry in item["items"]
                )
                chunks.append("</ul>")
            chunks.append("</dd>")
        chunks.append("</dl>")
    for issue in display["issues"]:
        chunks.append(
            f"<p><strong>{html.escape(issue['label'])}</strong>: "
            f"{html.escape(issue['value'])}</p>"
        )
    for reason in display["unknown_reasons"]:
        chunks.append(f"<p><strong>Unknown reason</strong>: {html.escape(reason)}</p>")
    chunks.append("</div>")
    return "".join(chunks)


def format_validation_summary(validation: Dict[str, Any]) -> str:
    findings = validation.get("key_findings", []) or []
    if not findings:
        return str(validation.get("state") or "not_recorded")
    codes = ", ".join(
        f"{item.get('category')}:{item.get('code')}" for item in findings
    )
    return f"{validation.get('state', 'not_recorded')} ({codes})"


def format_audit_summary(audit: Dict[str, Any]) -> str:
    if not audit:
        return "No post-finalize audit has been run for the current route"
    grades = audit.get("buyability_grades") or {}
    summary = (
        f"status={audit.get('status', 'unknown')}; "
        f"buyability B0/B1/B2={grades.get('B0', 0)}/{grades.get('B1', 0)}/{grades.get('B2', 0)}; "
        f"safety S1 substances={audit.get('s1_substances', 0)}; "
        f"process-ready/incomplete steps={audit.get('conditions_ready_steps', 0)}/"
        f"{audit.get('conditions_incomplete_steps', 0)}"
    )
    reasons = audit.get("incompleteness_reasons", []) or []
    if reasons:
        summary += "; incomplete because: " + ", ".join(str(item) for item in reasons)
    return summary


def _review_context_text(review: Dict[str, Any]) -> str:
    if review.get("kind") == "node":
        location = str(review.get("node_id") or review.get("smiles") or "reviewed node")
        text = f"{review.get('review_mode_label') or 'Node review'} at {location}"
        if review.get("reason"):
            text += f". Reason: {review['reason']}"
        if review.get("instruction") and review.get("instruction") != review.get("reason"):
            text += f". Direction: {review['instruction']}"
        return text
    text = str(review.get("smiles") or "Reviewed terminal")
    if review.get("reason"):
        text += f". Reason: {review['reason']}"
    return text


def _terminal_basis_text(basis: Dict[str, Any]) -> str:
    text = str(basis.get("source") or "Terminal decision recorded")
    if basis.get("reason"):
        text += f". Reason: {basis['reason']}"
    return text


def _append_text_scouting(lines: List[str], rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    lines.extend(["", "Scouting review summary:"])
    for row in rows:
        lines.append(
            f"  Node {row['node_id']} ({row['round_count']} round): "
            f"{row['review_summary']}"
        )
        for item in row["shortlist"]:
            lines.append(f"    shortlist: {item['summary']}")
            if item.get("frontier"):
                lines.append(
                    "      FRONTIER unprecedented_hypothesis: "
                    + (item["frontier"].get("basis") or "basis not supplied")
                )
        for item in row["deferred_review_seeds"]:
            lines.append(f"    review_node alternative: {item['summary']}")
            if item.get("frontier"):
                lines.append(
                    "      FRONTIER unprecedented_hypothesis: "
                    + (item["frontier"].get("basis") or "basis not supplied")
                )
        lines.append(
            f"    adoption summary: {row['adopted_count']} committed direction(s)"
            if row["adopted_count"]
            else "    reviewed; no scouting direction adopted"
        )


def generate_forward_report(
    tree: RetrosynthesisTree,
    scouting_view: Optional[Dict[str, Any]] = None,
    *,
    report_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
    route_plan: Optional[Dict[str, Any]] = None,
    experimental_outcomes: Optional[Iterable[Any]] = None,
    variant_info: Optional[Dict[str, Any]] = None,
    decision_history: Optional[Iterable[Any]] = None,
) -> str:
    view = report_view or build_route_report_view(
        tree,
        route_plan=route_plan,
        scouting_view=scouting_view,
        post_route_audit=post_route_audit,
        experimental_outcomes=experimental_outcomes,
        variant_info=variant_info,
        decision_history=decision_history,
    )
    lines = [
        f"Forward Synthesis Planning Report: {view['target_name']}",
        "=" * 54,
        f"Target: {view['target']}",
        f"Route scale: {view['step_count']} steps; "
        f"{view['starting_material_count']} planned starting materials",
    ]
    plan = view.get("route_plan") or {}
    if plan.get("route_thesis"):
        lines.append(f"Route thesis: {plan['route_thesis']}")
    if plan.get("route_mode"):
        lines.append(f"Route mode: {plan['route_mode']}")
    if view.get("route_summary"):
        lines.append(f"Route summary: {view['route_summary']}")
    if view.get("route_reviews"):
        lines.append("Route review context:")
        for review in view["route_reviews"]:
            lines.append(
                f"  - {review.get('title') or 'Route review'}: "
                + _review_context_text(review)
            )
    gate_counts = view.get("validation_gate_counts") or {}
    if gate_counts:
        lines.append(
            "Validation gate: "
            + ", ".join(f"{state}={count}" for state, count in gate_counts.items())
        )
    lines.append(
        "Experimental outcomes: "
        + (
            ", ".join(
                f"{status}={count}"
                for status, count in (view.get("experimental_outcome_counts") or {}).items()
            )
            if view.get("experimental_outcome_counts")
            else "none; planning commits do not imply experimental success"
        )
    )
    audit = view.get("post_route_audit") or {}
    lines.append("Post-route audit summary: " + format_audit_summary(audit))
    if audit:
        lines.append("Full audit: post_route_audit.md / post_route_audit.json")

    _append_text_scouting(lines, view.get("scouting_reviews", []))
    planning_stops = view.get("planning_stops", [])
    if planning_stops:
        lines.extend([
            "",
            f"Unresolved Planning Stops ({len(planning_stops)}):",
            "  These molecules were skipped to continue planning; they are not accepted starting materials and block finalization.",
        ])
        for item in planning_stops:
            lines.append(f"  - {item['smiles']}")
            if item.get("terminal_basis"):
                lines.append(
                    "    Planning-stop basis: "
                    + _terminal_basis_text(item["terminal_basis"])
                )
    lines.extend(["", f"Planned starting materials ({view['starting_material_count']}):"])
    for item in view["starting_materials"]:
        lines.append(f"  - {item['smiles']}")
        if item.get("terminal_basis"):
            lines.append(
                "    Terminal basis: "
                + _terminal_basis_text(item["terminal_basis"])
            )

    lines.extend(["", f"Forward synthesis steps ({view['step_count']}):", "-" * 54])
    for step in view["steps"]:
        lines.extend([
            "",
            f"Step {step['forward_step']} [{step['step_id']}]: {step['reaction_name']}",
            f"  Reaction: {step['reaction_smiles']}",
            "  Precursors: " + " + ".join(item["smiles"] for item in step["reactants"]),
        ])
        if step["reagents"]:
            lines.append("  Reagents for this step: " + " + ".join(step["reagents"]))
        lines.extend([
            f"  Product: {step['product']['smiles']}",
            "  Validation: " + format_validation_summary(step["validation"]),
        ])
        if step["validation_override"]:
            lines.append(
                "  Validation override: "
                + str(step["validation_override"].get("reason") or "declared")
            )
        if step["reasoning"]:
            lines.append(f"  Chemistry rationale: {step['reasoning']}")
        lines.append(f"  Risk assessment: {step['risk_assessment']}")
        lines.append("  Process conditions:")
        lines.extend(_text_condition_lines(step["conditions"]))
        for outcome in step["experimental_outcomes"]:
            outcome_text = str(outcome.get("outcome") or "unknown")
            if outcome.get("yield_percent") is not None:
                outcome_text += f", yield={outcome['yield_percent']}%"
            lines.append(f"  Experimental outcome: {outcome_text}")
        for rejected in step["rejected_alternatives"]:
            lines.append(f"  Rejected alternative: {rejected}")
    lines.extend([
        "",
        "Complete planning and audit records: session.json / tree.json.",
    ])
    return "\n".join(lines)


def _asset_path(path: Optional[str], report_dir: Path) -> str:
    if not path:
        return ""
    try:
        return os.path.relpath(path, report_dir).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def _md_scouting(lines: List[str], rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    lines.append("## Scouting Review Summary\n")
    for row in rows:
        lines.append(f"### Node {row['node_id']}\n")
        lines.append(f"- **Review**: {row['review_summary']}")
        lines.append(f"- **Rounds**: {row['round_count']}")
        for item in row["shortlist"]:
            lines.append(f"- **Shortlist**: {item['summary']}")
            if item.get("frontier"):
                lines.append(
                    "- **FRONTIER `unprecedented_hypothesis`**: "
                    + (item["frontier"].get("basis") or "basis not supplied")
                )
        for item in row["deferred_review_seeds"]:
            lines.append(f"- **`review_node` alternative**: {item['summary']}")
            if item.get("frontier"):
                lines.append(
                    "- **FRONTIER `unprecedented_hypothesis`**: "
                    + (item["frontier"].get("basis") or "basis not supplied")
                )
        lines.append(
            f"- **Adoption summary**: {row['adopted_count']} direction(s) entered the committed route"
            if row["adopted_count"]
            else "- **Adoption summary**: reviewed; no direction adopted"
        )
        lines.append("")


def write_markdown_report(
    tree: RetrosynthesisTree,
    mol_images: Dict[str, str],
    rxn_images: Dict[str, str],
    tree_image: Optional[str],
    text_report: str,
    output_path: str,
    mol_name: str,
    scouting_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
    report_view: Optional[Dict[str, Any]] = None,
    variant_info: Optional[Dict[str, Any]] = None,
    decision_history: Optional[Iterable[Any]] = None,
) -> None:
    del text_report
    view = report_view or build_route_report_view(
        tree,
        scouting_view=scouting_view,
        post_route_audit=post_route_audit,
        variant_info=variant_info,
        decision_history=decision_history,
    )
    report_dir = Path(output_path).resolve().parent
    lines = [
        f"# {mol_name} - Forward Synthesis Planning Report\n",
        f"**Target SMILES**: `{view['target']}`\n",
        f"**Route scale**: {view['step_count']} steps; "
        f"{view['starting_material_count']} planned starting materials\n",
    ]
    target = tree.get_molecule_by_smiles(tree.target)
    if target and target.node_id in mol_images:
        lines.append(
            f"![Target structure]({_asset_path(mol_images[target.node_id], report_dir)})\n"
        )
    plan = view.get("route_plan") or {}
    lines.append("## Route Conclusion\n")
    if plan.get("route_thesis"):
        lines.append(f"- **Route thesis**: {plan['route_thesis']}")
    if plan.get("route_mode"):
        lines.append(f"- **Route mode**: {plan['route_mode']}")
    if view.get("route_summary"):
        lines.append(f"- **Route summary**: {view['route_summary']}")
    if view.get("route_reviews"):
        lines.append("- **Route review context**:")
        for review in view["route_reviews"]:
            lines.append(
                f"  - **{review.get('title') or 'Route review'}**: "
                + _review_context_text(review)
            )
    gate_counts = view.get("validation_gate_counts") or {}
    if gate_counts:
        lines.append(
            "- **Validation gate**: "
            + ", ".join(f"`{state}`={count}" for state, count in gate_counts.items())
        )
    lines.append(
        "- **Experimental outcomes**: "
        + (
            ", ".join(
                f"`{status}`={count}"
                for status, count in (view.get("experimental_outcome_counts") or {}).items()
            )
            if view.get("experimental_outcome_counts")
            else "none; planning commits do not imply experimental success"
        )
    )
    audit = view.get("post_route_audit") or {}
    lines.append(f"- **Post-route audit summary**: {format_audit_summary(audit)}")
    if audit:
        lines.append(
            "- **Full audit**: [Markdown](post_route_audit.md) | "
            "[JSON](post_route_audit.json)"
        )
    lines.append("")
    if tree_image:
        lines.extend([
            "## Route Overview\n",
            f"![Route tree]({_asset_path(tree_image, report_dir)})\n",
        ])
    _md_scouting(lines, view.get("scouting_reviews", []))

    planning_stops = view.get("planning_stops", [])
    if planning_stops:
        lines.append(f"## Unresolved Planning Stops ({len(planning_stops)})\n")
        lines.append(
            "> These molecules were skipped to continue planning. They are not accepted "
            "starting materials and block finalization.\n"
        )
        for item in planning_stops:
            basis = item.get("terminal_basis") or {}
            lines.append(
                f"- `{item['smiles']}` — " + _terminal_basis_text(basis)
            )
        lines.append("")

    lines.append(f"## Planned Starting Materials ({view['starting_material_count']})\n")
    lines.extend(["| Structure | SMILES |", "|---|---|"])
    for item in view["starting_materials"]:
        image_path = mol_images.get(item["node_id"])
        image = (
            f"![{item['node_id']}]({_asset_path(image_path, report_dir)})"
            if image_path else "-"
        )
        lines.append(f"| {image} | `{item['smiles']}` |")
    lines.append("")
    for item in view["starting_materials"]:
        if item.get("terminal_basis"):
            lines.append(
                f"- **Terminal basis for `{item['node_id']}`**: "
                + _terminal_basis_text(item["terminal_basis"])
            )
    if any(item.get("terminal_basis") for item in view["starting_materials"]):
        lines.append("")

    lines.append(f"## Forward Synthesis Steps ({view['step_count']})\n")
    for step in view["steps"]:
        lines.append(
            f"### Step {step['forward_step']}: {step['reaction_name']} "
            f"(`{step['step_id']}`)\n"
        )
        image_path = rxn_images.get(step["step_id"])
        if image_path:
            lines.append(
                f"![Step {step['forward_step']}]({_asset_path(image_path, report_dir)})\n"
            )
        lines.extend([
            f"- **Reaction**: `{step['reaction_smiles']}`",
            "- **Precursors**: " + " + ".join(f"`{x['smiles']}`" for x in step["reactants"]),
        ])
        if step["reagents"]:
            lines.append("- **Reagents for this step**: " + " + ".join(f"`{x}`" for x in step["reagents"]))
        lines.extend([
            f"- **Product**: `{step['product']['smiles']}`",
            f"- **Validation**: {format_validation_summary(step['validation'])}",
        ])
        if step["validation_override"]:
            lines.append(
                "- **Validation override**: "
                + str(step["validation_override"].get("reason") or "declared")
            )
        if step["reasoning"]:
            lines.append(f"- **Chemistry rationale**: {step['reasoning']}")
        lines.append(f"- **Risk assessment**: {step['risk_assessment']}")
        lines.extend(_markdown_condition_lines(step["conditions"]))
        for outcome in step["experimental_outcomes"]:
            lines.append(f"- **Experimental outcome**: `{outcome.get('outcome', 'unknown')}`")
        for rejected in step["rejected_alternatives"]:
            lines.append(f"- **Rejected alternative**: {rejected}")
        lines.append("")
    lines.extend([
        "## Raw Data\n",
        "Complete planning and audit records are retained in [`session.json`](session.json) "
        "and [`tree.json`](tree.json).\n",
    ])
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")


_HTML_STYLE = """
*{box-sizing:border-box}body{margin:0;background:#f4f5f3;color:#20231f;font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1080px;margin:auto;padding:24px}header,section{background:#fff;border:1px solid #d9ddd7;border-radius:8px;padding:20px;margin-bottom:16px}
h1{font-size:24px;margin:0 0 8px}h2{font-size:18px;margin:0 0 14px;border-bottom:1px solid #e4e7e2;padding-bottom:8px}h3{font-size:16px;margin:0 0 10px}
.muted{color:#60675f}.status{border-left:4px solid #477a54;background:#f3f8f3;padding:10px 12px}.warning{border-left:4px solid #ad791d;background:#fff9eb;padding:10px 12px}
.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin-top:14px}.metric{border:1px solid #e1e4df;padding:10px;border-radius:6px}
.target{max-width:320px;max-height:220px;display:block;margin:14px auto}.tree,.reaction{max-width:100%;display:block;margin:10px auto}.route-shell{border:1px solid #d8ddd7;border-radius:6px;background:#f7f8f6;overflow:hidden}.route-toolbar{min-height:46px;display:flex;align-items:center;justify-content:space-between;gap:10px;padding:7px 8px;border-bottom:1px solid #d8ddd7;background:#fff}.route-tools{display:flex;align-items:center;gap:6px;flex-wrap:wrap}.route-toolbar button{min-width:36px;height:32px;border:1px solid #cbd1ca;border-radius:6px;background:#fff;color:#292d28;font:600 13px/1 system-ui,-apple-system,"Segoe UI",sans-serif;letter-spacing:0;padding:0 10px;cursor:pointer}.route-toolbar button:hover{background:#f0f3ef}.route-toolbar button:focus-visible{outline:2px solid #477a54;outline-offset:2px}.route-toolbar button[aria-pressed=true]{background:#eaf2ec;border-color:#78947f}.route-zoom{min-width:48px;color:#596058;font-variant-numeric:tabular-nums;text-align:right}.route-viewport{height:clamp(420px,72vh,760px);overflow:auto;background:#fff;position:relative;cursor:grab;touch-action:none}.route-viewport.is-panning{cursor:grabbing}.route-diagram{display:block;max-width:none;height:auto;margin:0 auto;user-select:none;-webkit-user-drag:none}.route-shell:fullscreen{height:100vh;border:0;border-radius:0;background:#fff}.route-shell:fullscreen .route-viewport{height:calc(100vh - 46px)}.materials{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px}.material{border:1px solid #e1e4df;padding:10px;border-radius:6px}.material img{max-width:100%;height:130px;object-fit:contain}
.step{border-top:1px solid #dde1db;padding:18px 0}.step:first-of-type{border-top:0}.facts{display:grid;grid-template-columns:140px 1fr;gap:6px 12px}.label{font-weight:650}.detail-panel{padding:10px 12px;border:1px solid;border-radius:6px}.rationale{background:#f4f7fa;border-color:#dfe8ee}.risk{background:#fff9eb;border-color:#efe1bd}.conditions{background:#f4f8f6;border-color:#dfe9e3}.conditions>p{margin:0 0 8px}.conditions h4{font-size:14px;margin:12px 0 6px}.conditions dl{display:grid;grid-template-columns:120px 1fr;gap:5px 10px;margin:0}.conditions dt{font-weight:650}.conditions dd{margin:0}.conditions ol,.conditions ul{margin:0;padding-left:22px}.review-context{background:#f4f7fa;border:1px solid #dfe8ee;border-left:4px solid #718899;border-radius:6px;padding:10px 12px;margin-top:12px}.review-context p{margin:4px 0}.terminal-note{background:#f6f7f5;border:1px solid #e0e3df;border-radius:6px;padding:8px 10px;margin-top:10px;font-size:13px}.frontier{background:#fff0ed;border-left:4px solid #b74432;padding:10px}.mono,code{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;word-break:break-word}a{color:#2f6340}.raw{font-size:13px;color:#596058}
@media(max-width:640px){main{padding:10px}.facts{grid-template-columns:1fr}.label{margin-top:6px}.route-toolbar{align-items:flex-start}.route-viewport{height:clamp(380px,68vh,620px)}}
"""


_ROUTE_VIEWER_SCRIPT = r"""
(() => {
  const clamp = (value, minimum, maximum) => Math.min(maximum, Math.max(minimum, value));
  document.querySelectorAll('[data-route-viewer]').forEach((viewer) => {
    const viewport = viewer.querySelector('[data-route-viewport]');
    const diagram = viewer.querySelector('[data-route-diagram]');
    const zoomOutput = viewer.querySelector('[data-route-zoom]');
    const fullscreenButton = viewer.querySelector('[data-route-action="fullscreen"]');
    let naturalWidth = 0;
    let naturalHeight = 0;
    let scale = 1;
    let mode = 'fit-width';
    let drag = null;

    const updateButtons = () => {
      viewer.querySelectorAll('[data-route-action="fit-width"],[data-route-action="fit-route"]').forEach((button) => {
        button.setAttribute('aria-pressed', String(button.dataset.routeAction === mode));
      });
      zoomOutput.textContent = `${Math.round(scale * 100)}%`;
    };

    const setScale = (nextScale, nextMode = 'custom') => {
      if (!naturalWidth || !naturalHeight) return;
      const previousWidth = naturalWidth * scale;
      const previousHeight = naturalHeight * scale;
      const centerX = (viewport.scrollLeft + (viewport.clientWidth / 2)) / Math.max(previousWidth, 1);
      const centerY = (viewport.scrollTop + (viewport.clientHeight / 2)) / Math.max(previousHeight, 1);
      scale = clamp(nextScale, 0.2, 4);
      mode = nextMode;
      diagram.style.width = `${Math.round(naturalWidth * scale)}px`;
      diagram.style.height = 'auto';
      viewport.scrollLeft = (centerX * naturalWidth * scale) - (viewport.clientWidth / 2);
      viewport.scrollTop = (centerY * naturalHeight * scale) - (viewport.clientHeight / 2);
      updateButtons();
    };

    const fitWidth = () => {
      const availableWidth = Math.max(1, viewport.clientWidth - 16);
      setScale(Math.min(1, availableWidth / naturalWidth), 'fit-width');
      viewport.scrollTop = 0;
    };

    const fitRoute = () => {
      const availableWidth = Math.max(1, viewport.clientWidth - 16);
      const availableHeight = Math.max(1, viewport.clientHeight - 16);
      setScale(Math.min(1, availableWidth / naturalWidth, availableHeight / naturalHeight), 'fit-route');
      viewport.scrollLeft = 0;
      viewport.scrollTop = 0;
    };

    const initialize = () => {
      naturalWidth = diagram.naturalWidth;
      naturalHeight = diagram.naturalHeight;
      fitWidth();
    };

    viewer.addEventListener('click', (event) => {
      const button = event.target.closest('[data-route-action]');
      if (!button) return;
      const action = button.dataset.routeAction;
      if (action === 'zoom-in') setScale(scale * 1.2);
      if (action === 'zoom-out') setScale(scale / 1.2);
      if (action === 'actual-size') setScale(1);
      if (action === 'fit-width') fitWidth();
      if (action === 'fit-route') fitRoute();
      if (action === 'fullscreen') {
        if (document.fullscreenElement === viewer) {
          document.exitFullscreen();
        } else if (viewer.requestFullscreen) {
          viewer.requestFullscreen().catch(() => {});
        }
      }
    });

    viewport.addEventListener('pointerdown', (event) => {
      if (event.button !== 0) return;
      drag = {
        x: event.clientX,
        y: event.clientY,
        left: viewport.scrollLeft,
        top: viewport.scrollTop,
      };
      viewport.classList.add('is-panning');
      viewport.setPointerCapture(event.pointerId);
      event.preventDefault();
    });
    viewport.addEventListener('pointermove', (event) => {
      if (!drag) return;
      viewport.scrollLeft = drag.left - (event.clientX - drag.x);
      viewport.scrollTop = drag.top - (event.clientY - drag.y);
    });
    const endDrag = (event) => {
      if (!drag) return;
      drag = null;
      viewport.classList.remove('is-panning');
      if (viewport.hasPointerCapture(event.pointerId)) viewport.releasePointerCapture(event.pointerId);
    };
    viewport.addEventListener('pointerup', endDrag);
    viewport.addEventListener('pointercancel', endDrag);
    diagram.addEventListener('dragstart', (event) => event.preventDefault());
    diagram.addEventListener('load', initialize);

    document.addEventListener('fullscreenchange', () => {
      const isFullscreen = document.fullscreenElement === viewer;
      fullscreenButton.textContent = isFullscreen ? 'Exit full screen' : 'Full screen';
      fullscreenButton.setAttribute('aria-label', fullscreenButton.textContent);
      if (mode === 'fit-width') fitWidth();
      if (mode === 'fit-route') fitRoute();
    });
    window.addEventListener('resize', () => {
      if (mode === 'fit-width') fitWidth();
      if (mode === 'fit-route') fitRoute();
    });
    if (diagram.complete) initialize();
  });
})();
"""


def _html_scouting(rows: List[Dict[str, Any]]) -> str:
    if not rows:
        return ""
    chunks = ["<section><h2>Scouting Review Summary</h2>"]
    for row in rows:
        chunks.append(
            f"<div class='step'><h3>Node {html.escape(row['node_id'])}</h3>"
            f"<p>{html.escape(row['review_summary'])}</p>"
            f"<p class='muted'>Rounds: {row['round_count']}; "
            f"committed adoption count: {row['adopted_count']}</p>"
        )
        for label, items in (
            ("Shortlist", row["shortlist"]),
            ("review_node alternative", row["deferred_review_seeds"]),
        ):
            for item in items:
                chunks.append(f"<p><strong>{label}</strong>: {html.escape(item['summary'])}</p>")
                if item.get("frontier"):
                    basis = item["frontier"].get("basis") or "basis not supplied"
                    chunks.append(
                        "<p class='frontier'><strong>FRONTIER "
                        "unprecedented_hypothesis</strong>: "
                        f"{html.escape(basis)}</p>"
                    )
        if not row["adopted_count"]:
            chunks.append("<p class='muted'>Reviewed; no scouting direction adopted.</p>")
        chunks.append("</div>")
    chunks.append("</section>")
    return "".join(chunks)


def generate_html_report(
    tree: RetrosynthesisTree,
    mol_name: str,
    tree_image_path: Optional[str] = None,
    scouting_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
    report_view: Optional[Dict[str, Any]] = None,
    variant_info: Optional[Dict[str, Any]] = None,
    decision_history: Optional[Iterable[Any]] = None,
    tree_svg_path: Optional[str] = None,
) -> str:
    view = report_view or build_route_report_view(
        tree,
        scouting_view=scouting_view,
        post_route_audit=post_route_audit,
        variant_info=variant_info,
        decision_history=decision_history,
    )
    target = tree.get_molecule_by_smiles(tree.target)
    target_src = f"images/{target.node_id}.png" if target else ""
    chunks = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>{html.escape(mol_name)} - Forward Synthesis Planning Report</title>",
        f"<style>{_HTML_STYLE}</style></head><body><main>",
        f"<header><h1>{html.escape(mol_name)} - Forward Synthesis Planning Report</h1>",
        f"<p class='mono'>{html.escape(view['target'])}</p>",
    ]
    if target_src:
        chunks.append(f"<img class='target' src='{target_src}' alt='Target structure'>")
    chunks.append("<div class='summary'>")
    chunks.append(f"<div class='metric'><strong>Steps</strong><br>{view['step_count']}</div>")
    chunks.append(
        f"<div class='metric'><strong>Planned starting materials</strong><br>{view['starting_material_count']}</div>"
    )
    gate_counts = view.get("validation_gate_counts") or {}
    chunks.append(
        "<div class='metric'><strong>Validation</strong><br>"
        + html.escape(", ".join(f"{k}={v}" for k, v in gate_counts.items()) or "not recorded")
        + "</div>"
    )
    outcomes = view.get("experimental_outcome_counts") or {}
    chunks.append(
        "<div class='metric'><strong>Experimental outcomes</strong><br>"
        + html.escape(", ".join(f"{k}={v}" for k, v in outcomes.items()) or "none recorded")
        + "</div></div></header>"
    )

    plan = view.get("route_plan") or {}
    chunks.append("<section><h2>Route Conclusion</h2>")
    if plan.get("route_thesis"):
        chunks.append(f"<p><strong>Route thesis</strong>: {html.escape(plan['route_thesis'])}</p>")
    if plan.get("route_mode"):
        chunks.append(f"<p><strong>Route mode</strong>: {html.escape(plan['route_mode'])}</p>")
    if view.get("route_summary"):
        chunks.append(f"<p><strong>Route summary</strong>: {html.escape(view['route_summary'])}</p>")
    if view.get("route_reviews"):
        chunks.append("<div class='review-context'><strong>Route review context</strong>")
        for review in view["route_reviews"]:
            chunks.append(
                f"<p><strong>{html.escape(str(review.get('title') or 'Route review'))}</strong>: "
                f"{html.escape(_review_context_text(review))}</p>"
            )
        chunks.append("</div>")
    audit = view.get("post_route_audit") or {}
    chunks.append(f"<p><strong>Post-route audit summary</strong>: {html.escape(format_audit_summary(audit))}</p>")
    if audit:
        chunks.append(
            "<p><a href='post_route_audit.md'>Full audit Markdown</a> | "
            "<a href='post_route_audit.json'>Machine-readable JSON</a></p>"
        )
    chunks.append("</section>")
    if tree_svg_path or tree_image_path:
        route_src = (
            "images/synthesis_tree.svg"
            if tree_svg_path
            else "images/synthesis_tree.png"
        )
        fallback = (
            " onerror=\"this.onerror=null;this.src='images/synthesis_tree.png'\""
            if tree_svg_path and tree_image_path
            else ""
        )
        chunks.append(
            "<section><h2>Route Overview</h2>"
            "<div class='route-shell' data-route-viewer>"
            "<div class='route-toolbar' role='toolbar' aria-label='Route overview controls'>"
            "<div class='route-tools'>"
            "<button type='button' data-route-action='zoom-out' title='Zoom out' aria-label='Zoom out'>-</button>"
            "<button type='button' data-route-action='zoom-in' title='Zoom in' aria-label='Zoom in'>+</button>"
            "<button type='button' data-route-action='actual-size' title='Actual size'>100%</button>"
            "<button type='button' data-route-action='fit-width' title='Fit width' aria-pressed='true'>Fit width</button>"
            "<button type='button' data-route-action='fit-route' title='Fit whole route' aria-pressed='false'>Fit route</button>"
            "</div>"
            "<div class='route-tools'>"
            "<output class='route-zoom' data-route-zoom aria-live='polite'>100%</output>"
            "<button type='button' data-route-action='fullscreen' title='Full screen' aria-label='Full screen'>Full screen</button>"
            "</div></div>"
            "<div class='route-viewport' data-route-viewport tabindex='0' aria-label='Scrollable route overview'>"
            f"<img class='route-diagram' data-route-diagram src='{route_src}'{fallback} "
            "alt='Target at top with reaction steps and precursors expanding downward'>"
            "</div></div></section>"
        )
    chunks.append(_html_scouting(view.get("scouting_reviews", [])))

    planning_stops = view.get("planning_stops", [])
    if planning_stops:
        chunks.append(
            f"<section class='warning'><h2>Unresolved Planning Stops ({len(planning_stops)})</h2>"
            "<p>These molecules were skipped to continue planning. They are not accepted "
            "starting materials and block finalization.</p><ul>"
        )
        for item in planning_stops:
            basis = item.get("terminal_basis") or {}
            chunks.append(
                f"<li><code>{html.escape(item['smiles'])}</code> — "
                + html.escape(_terminal_basis_text(basis))
                + "</li>"
            )
        chunks.append("</ul></section>")

    chunks.append("<section><h2>Planned Starting Materials</h2><div class='materials'>")
    for item in view["starting_materials"]:
        chunks.append(
            "<div class='material'>"
            f"<img src='images/{html.escape(item['node_id'])}.png' alt='Starting material'>"
            f"<code>{html.escape(item['smiles'])}</code>"
        )
        if item.get("terminal_basis"):
            chunks.append(
                "<div class='terminal-note'><strong>Terminal basis</strong>: "
                + html.escape(_terminal_basis_text(item["terminal_basis"]))
                + "</div>"
            )
        chunks.append("</div>")
    chunks.append("</div></section>")

    chunks.append("<section><h2>Forward Synthesis Steps</h2>")
    for step in view["steps"]:
        chunks.append(
            f"<article class='step'><h3>Step {step['forward_step']}: "
            f"{html.escape(step['reaction_name'])} "
            f"<span class='muted'>({html.escape(step['step_id'])})</span></h3>"
            f"<img class='reaction' src='images/{html.escape(step['step_id'])}_reaction.png' "
            f"alt='Step {step['forward_step']} reaction'>"
            "<div class='facts'>"
            f"<div class='label'>Precursors</div><div>{html.escape(' + '.join(x['smiles'] for x in step['reactants']))}</div>"
        )
        if step["reagents"]:
            chunks.append(
                f"<div class='label'>Reagents for this step</div><div>{html.escape(' + '.join(step['reagents']))}</div>"
            )
        chunks.extend([
            f"<div class='label'>Product</div><div>{html.escape(step['product']['smiles'])}</div>",
            f"<div class='label'>Validation</div><div>{html.escape(format_validation_summary(step['validation']))}</div>",
        ])
        if step["validation_override"]:
            chunks.append(
                "<div class='label'>Validation override</div><div class='warning'>"
                + html.escape(str(step["validation_override"].get("reason") or "declared"))
                + "</div>"
            )
        if step["reasoning"]:
            chunks.append(
                f"<div class='label'>Chemistry rationale</div><div class='rationale detail-panel'>{html.escape(step['reasoning'])}</div>"
            )
        chunks.extend([
            f"<div class='label'>Risk assessment</div><div class='risk detail-panel'>{html.escape(step['risk_assessment'])}</div>",
            f"<div class='label'>Process conditions</div>{_html_condition_block(step['conditions'])}",
        ])
        for outcome in step["experimental_outcomes"]:
            chunks.append(
                f"<div class='label'>Experimental outcome</div><div>{html.escape(str(outcome.get('outcome') or 'unknown'))}</div>"
            )
        for rejected in step["rejected_alternatives"]:
            chunks.append(
                f"<div class='label'>Rejected alternative</div><div>{html.escape(rejected)}</div>"
            )
        chunks.append("</div></article>")
    chunks.append("</section>")
    chunks.append(
        "<section class='raw'><h2>Raw Data</h2>"
        "<p>Complete planning and audit records are retained in "
        "<a href='session.json'>session.json</a> and "
        "<a href='tree.json'>tree.json</a>.</p></section></main>"
    )
    if tree_svg_path or tree_image_path:
        chunks.append(f"<script>{_ROUTE_VIEWER_SCRIPT}</script>")
    chunks.append("</body></html>")
    return "".join(chunks)
