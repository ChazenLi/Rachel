"""Independent post-finalize buyability and safety audit for Rachel routes."""

from __future__ import annotations

import hashlib
import html
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from rdkit import Chem

from Rachel.chem_tools._rdkit_utils import load_template, parse_mol, smarts_match
from Rachel.chem_tools.fg_warnings import check_fg_conflicts
from Rachel.display_language import dangerous_combo_text
from Rachel.main.report_projection import forward_reactions
from Rachel.main.retro_report import get_terminal_list
from Rachel.main.risk_protocol import (
    PROCESS_CONDITIONS_INSUFFICIENT,
    PROCESS_CONDITIONS_PRELIMINARY_REVIEW,
    assess_process_conditions,
    format_process_condition_gaps,
)
from Rachel.tools.pubchem_terminal_audit import (
    PUBCHEM_BASE,
    PubChemClient,
    audit_record,
    summarize as summarize_buyability,
)


POST_ROUTE_AUDIT_SCHEMA = "rachel.post_route_audit.v2"
DEFAULT_SCOPES = ("buyability", "safety")


def route_digest(tree: Any) -> str:
    payload = json.dumps(
        tree.to_dict(),
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return f"sha256:{hashlib.sha256(payload).hexdigest()}"


def audit_matches_tree(payload: Dict[str, Any], tree: Any) -> bool:
    return bool(payload) and payload.get("route_digest") == route_digest(tree)


def _unique_strings(values: Iterable[Any]) -> List[str]:
    result: List[str] = []
    seen = set()
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            result.append(text)
    return result


def _inchikey(smiles: str) -> str:
    mol = parse_mol(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


def _identity_match_level(
    input_smiles: str,
    query_smiles: str,
    properties: Dict[str, Any],
    normalization_notes: Optional[List[str]] = None,
) -> str:
    if normalization_notes or input_smiles != query_smiles:
        normalized = True
    else:
        normalized = False
    input_key = _inchikey(query_smiles)
    returned_key = str(properties.get("InChIKey", "") or "")
    if input_key and returned_key and input_key == returned_key:
        return "normalized_parent" if normalized else "exact_full"
    if input_key and returned_key and input_key.split("-")[0] == returned_key.split("-")[0]:
        return "normalized_parent" if normalized else "connectivity_only"
    return "unresolved"


def _buyability_section(
    tree: Any,
    client: PubChemClient,
    *,
    include_vendors: bool,
    query_reagents: bool,
) -> Dict[str, Any]:
    results = [
        audit_record(
            record,
            client,
            include_vendors=include_vendors,
            query_reagents=query_reagents,
        )
        for record in get_terminal_list(tree)
    ]
    for item in results:
        pubchem = item.get("pubchem") or {}
        metrics = item.get("pubchem_metrics") or {}
        local = item.get("local") or {}
        vendor = pubchem.get("vendor_evidence") or {}
        catalog_status = str(vendor.get("vendor_catalog_status", "not_queried"))
        if catalog_status == "hit":
            grade = "B2"
        elif metrics.get("cid_hit"):
            grade = "B1"
        else:
            grade = "B0"
        item["evidence_grade"] = grade
        item["identity_match_level"] = _identity_match_level(
            str(item.get("smiles", "")),
            str(pubchem.get("query_smiles", "")),
            pubchem.get("properties") or {},
            local.get("normalization_notes") or [],
        ) if metrics.get("cid_hit") else "unresolved"
        item["vendor_catalog"] = {
            "status": catalog_status,
            "vendor_count": vendor.get("vendor_count"),
            "listing_count": vendor.get("vendor_listing_count"),
            "vendor_names": list(vendor.get("vendor_names") or []),
            "rows": list(vendor.get("vendor_rows") or []),
            "rows_truncated": bool(vendor.get("vendor_rows_truncated")),
            "retrieved_at": vendor.get("retrieved_at", ""),
            "from_cache": bool(vendor.get("from_cache")),
        }
    summary = summarize_buyability(results)
    summary["by_evidence_grade"] = {
        grade: sum(1 for item in results if item.get("evidence_grade") == grade)
        for grade in ("B0", "B1", "B2")
    }
    return {
        "policy": "pubchem_identity_and_vendor_catalog_only",
        "summary": summary,
        "terminals": results,
        "limitations": [
            "B1 proves a PubChem identity record, not commercial availability.",
            "B2 proves depositor/vendor catalog rows, not live stock, price, lead time, region, purity, pack size, or procurement.",
            "No hit means no evidence in the checked PubChem sources; it does not prove that a material cannot be purchased.",
        ],
    }


def _material_records(tree: Any) -> List[Dict[str, Any]]:
    records: Dict[str, Dict[str, Any]] = {}

    def add(smiles: str, role: str, *, node_id: str = "", step_id: str = "") -> None:
        raw = str(smiles or "").strip()
        if not raw:
            return
        mol = parse_mol(raw)
        key = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) if mol else raw
        item = records.setdefault(
            key,
            {
                "entity_id": f"material_{len(records) + 1}",
                "smiles": raw,
                "canonical_smiles": key if mol else "",
                "roles": [],
                "node_ids": [],
                "step_ids": [],
            },
        )
        item["roles"] = _unique_strings(list(item["roles"]) + [role])
        item["node_ids"] = _unique_strings(list(item["node_ids"]) + [node_id])
        item["step_ids"] = _unique_strings(list(item["step_ids"]) + [step_id])

    for node in tree.molecule_nodes.values():
        add(node.smiles, str(node.role or "molecule"), node_id=node.node_id)
    for reaction in tree.reaction_nodes:
        for reagent in reaction.reagents:
            add(reagent, "reagent", step_id=reaction.step_id)
    return list(records.values())


def _structural_alerts(smiles: str, knowledge_profile: Any) -> List[Dict[str, Any]]:
    mol = parse_mol(smiles)
    if mol is None:
        return []
    alerts = load_template("structural_alerts.json", knowledge_profile=knowledge_profile)
    results: List[Dict[str, Any]] = []
    for alert_id, smarts in alerts.items():
        if not isinstance(smarts, str):
            continue
        matches = smarts_match(mol, smarts)
        if not matches:
            continue
        results.append({
            "alert_id": alert_id,
            "match_count": len(matches),
            "atom_matches": [list(match) for match in matches[:5]],
            "evidence_grade": "S0",
            "basis": "local_structural_screen",
            "interpretation": (
                "A low-specificity structural screening alert, not a toxicology "
                "classification or proof of hazard for this material."
            ),
            "knowledge_ref": knowledge_profile.source("chem.structural_alerts", alert_id),
        })
    return results


def _material_safety_record(
    record: Dict[str, Any],
    client: PubChemClient,
    knowledge_profile: Any,
    ghs_by_cid: Dict[int, Dict[str, Any]],
) -> Dict[str, Any]:
    result = dict(record)
    smiles = str(record.get("canonical_smiles") or record.get("smiles") or "")
    result["structural_screening"] = _structural_alerts(smiles, knowledge_profile)
    if parse_mol(smiles) is None:
        result.update({
            "identity_status": "invalid_or_non_smiles",
            "identity_match_level": "unresolved",
            "cid": None,
            "ghs": {"state": "not_queried", "classifications": []},
            "status": "incomplete",
            "evidence_grade": "S0",
        })
        return result

    lookup = client.cid_lookup(smiles)
    result["identity_status"] = lookup.get("status", "not_queried")
    result["cid_candidates"] = list(lookup.get("cids") or [])[:10]
    if lookup.get("status") != "hit" or not lookup.get("cids"):
        result.update({
            "identity_match_level": "unresolved",
            "cid": None,
            "ghs": {
                "state": "not_queried" if lookup.get("status") == "not_queried" else "identity_unresolved",
                "classifications": [],
            },
            "status": "incomplete",
            "evidence_grade": "S0",
        })
        if lookup.get("error"):
            result["identity_error"] = lookup.get("error")
        return result

    cid = int(lookup["cids"][0])
    properties = client.properties(cid)
    result["cid"] = cid
    result["identity"] = properties
    result["identity_match_level"] = _identity_match_level(
        smiles, smiles, properties, []
    )
    if cid not in ghs_by_cid:
        ghs_by_cid[cid] = client.ghs_summary(cid)
    result["ghs"] = ghs_by_cid[cid]
    if result["ghs"].get("state") == "classification_data_found":
        result["status"] = "flags_found"
        result["evidence_grade"] = "S1"
    elif result["ghs"].get("state") in {"not_queried", "lookup_error"}:
        result["status"] = "incomplete"
        result["evidence_grade"] = "S0"
    else:
        result["status"] = "no_classification_found_in_checked_sources"
        result["evidence_grade"] = "S0"
    return result


def _step_records(tree: Any, knowledge_profile: Any) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for forward_step, reaction in enumerate(forward_reactions(tree), 1):
        reactants = [
            tree.molecule_nodes[node_id].smiles
            for node_id in reaction.reactant_nodes
            if node_id in tree.molecule_nodes
        ]
        product = tree.molecule_nodes.get(reaction.product_node)
        components = reactants + list(reaction.reagents)
        combined = ".".join(value for value in components if parse_mol(value) is not None)
        screening = (
            check_fg_conflicts(combined, knowledge_profile=knowledge_profile)
            if combined
            else {"ok": False, "dangerous_combos": []}
        )
        dangerous = list(screening.get("dangerous_combos") or [])
        declaration = ""
        declared_conditions: Dict[str, Any] = {}
        if reaction.llm_decision is not None:
            declaration = str(reaction.llm_decision.risk_assessment or "").strip()
            conditions = getattr(reaction.llm_decision, "process_conditions", {})
            if isinstance(conditions, dict):
                declared_conditions = conditions
        condition_coverage = assess_process_conditions(declared_conditions)
        alerts = [
            {
                "category": "local_dangerous_combination_screen",
                "severity": "unknown",
                "basis": "rule",
                "rationale": str(item.get("warning", "") or ""),
                "groups": list(item.get("groups") or []),
                "evidence_grade": "S0",
                "required_review": "Confirm against exact reagents, conditions, SDS, literature, and laboratory SOP.",
            }
            for item in dangerous
        ]
        records.append({
            "forward_step": forward_step,
            "step_id": reaction.step_id,
            "reaction_name": reaction.reaction_type,
            "reactants": reactants,
            "reagents": list(reaction.reagents),
            "product": product.smiles if product else "",
            "declared_risk_assessment": declaration,
            "risk_assessment_status": "declared" if declaration else "not_assessed",
            "risk_assessment_evidence_type": "llm_authored_planning_judgment",
            "risk_assessment_semantic_validation": "not_performed",
            "alerts": alerts,
            "declared_conditions": declared_conditions,
            "process_condition_coverage": condition_coverage,
            "process_assessment_status": condition_coverage["process_assessment_status"],
            "evidence_grade": "S0",
        })
    return records


def _safety_section(tree: Any, client: PubChemClient, knowledge_profile: Any) -> Dict[str, Any]:
    ghs_by_cid: Dict[int, Dict[str, Any]] = {}
    substances = [
        _material_safety_record(record, client, knowledge_profile, ghs_by_cid)
        for record in _material_records(tree)
    ]
    steps = _step_records(tree, knowledge_profile)
    return {
        "summary": {
            "substance_count": len(substances),
            "substances_with_source_attributed_ghs": sum(
                1 for item in substances if item.get("evidence_grade") == "S1"
            ),
            "substances_with_structural_screening_alerts": sum(
                1 for item in substances if item.get("structural_screening")
            ),
            "step_count": len(steps),
            "steps_with_explicit_risk_assessment": sum(
                1 for item in steps if item.get("risk_assessment_status") == "declared"
            ),
            "steps_missing_explicit_risk_assessment": sum(
                1 for item in steps if item.get("risk_assessment_status") == "not_assessed"
            ),
            "steps_with_process_conditions_for_preliminary_review": sum(
                1
                for item in steps
                if item.get("process_assessment_status")
                == PROCESS_CONDITIONS_PRELIMINARY_REVIEW
            ),
            "steps_with_incomplete_process_conditions": sum(
                1
                for item in steps
                if item.get("process_assessment_status")
                == PROCESS_CONDITIONS_INSUFFICIENT
            ),
            "steps_with_local_screening_alerts": sum(
                1 for item in steps if item.get("alerts")
            ),
        },
        "substances": substances,
        "steps": steps,
        "limitations": [
            "GHS classifications are source-specific substance evidence and may differ across authorities or depositors.",
            "Supplier- and product-specific SDS records were not retrieved in this version.",
            "Structured condition coverage checks declared field presence and shape only; it does not validate chemical correctness or establish process safety.",
            "Even complete declared condition coverage supports only preliminary review and does not replace SDS, SOP, calorimetry, exposure, scale-up, or equipment assessment.",
            "S0 structural and combination alerts are screening hypotheses, not process hazard conclusions.",
            "Absence of a checked-source classification or alert does not mean safe.",
        ],
    }


def run_post_route_audit(
    tree: Any,
    knowledge_profile: Any,
    *,
    cache_dir: Optional[Path] = None,
    scopes: Sequence[str] = DEFAULT_SCOPES,
    offline: bool = False,
    include_vendors: bool = True,
    query_reagents: bool = False,
    timeout: float = 20.0,
    pause_seconds: float = 0.25,
) -> Dict[str, Any]:
    selected_scopes = _unique_strings(scopes)
    invalid_scopes = [scope for scope in selected_scopes if scope not in DEFAULT_SCOPES]
    if invalid_scopes:
        raise ValueError(f"unsupported audit scopes: {', '.join(invalid_scopes)}")
    generated_at = datetime.now(timezone.utc).isoformat()
    digest = route_digest(tree)
    client = PubChemClient(
        cache_dir=cache_dir,
        timeout=timeout,
        pause_seconds=max(0.2, pause_seconds),
        offline=offline,
    )
    payload: Dict[str, Any] = {
        "schema": POST_ROUTE_AUDIT_SCHEMA,
        "audit_id": f"audit_{generated_at.replace(':', '').replace('-', '')}_{digest[-8:]}",
        "generated_at": generated_at,
        "route_digest": digest,
        "scopes": selected_scopes,
        "offline": bool(offline),
        "status": "complete",
        "source_manifest": [
            {
                "source": "PubChem PUG REST",
                "url": f"{PUBCHEM_BASE}/pug",
                "purpose": "CID and identity properties",
            },
            {
                "source": "PubChem PUG View",
                "url": f"{PUBCHEM_BASE}/pug_view",
                "purpose": "vendor catalog rows and source-attributed GHS classifications",
            },
        ],
        "limitations": [
            "This audit does not change route chemistry, terminal decisions, validation gates, or route ranking.",
            "No result in this payload is a generic safe/unsafe verdict.",
        ],
    }
    incomplete_reasons: List[str] = []
    if "buyability" in selected_scopes:
        payload["buyability"] = _buyability_section(
            tree,
            client,
            include_vendors=include_vendors and not offline,
            query_reagents=query_reagents,
        )
        for item in payload["buyability"].get("terminals", []):
            metrics = item.get("pubchem_metrics") or {}
            catalog = item.get("vendor_catalog") or {}
            if metrics.get("cid_status") == "lookup_error":
                incomplete_reasons.append("pubchem_cid_lookup_error")
            if metrics.get("cid_hit") and item.get("identity_match_level") == "unresolved":
                incomplete_reasons.append("pubchem_identity_unresolved")
            if include_vendors and catalog.get("status") in {"lookup_error", "aggregate_flag_only"}:
                incomplete_reasons.append("pubchem_vendor_catalog_unresolved")
    if "safety" in selected_scopes:
        payload["safety"] = _safety_section(tree, client, knowledge_profile)
        safety_summary = payload["safety"]["summary"]
        if safety_summary.get("steps_missing_explicit_risk_assessment"):
            incomplete_reasons.append("steps_missing_explicit_risk_assessment")
        if any(
            item.get("status") == "incomplete"
            for item in payload["safety"].get("substances", [])
        ):
            incomplete_reasons.append("substance_safety_lookup_incomplete")
        if any(
            item.get("process_assessment_status")
            == PROCESS_CONDITIONS_INSUFFICIENT
            for item in payload["safety"].get("steps", [])
        ):
            incomplete_reasons.append("process_conditions_insufficient")
    if offline:
        incomplete_reasons.append("offline_mode")
    payload["incompleteness_reasons"] = _unique_strings(incomplete_reasons)
    if payload["incompleteness_reasons"]:
        payload["status"] = "incomplete"
    return payload


def _md_cell(value: Any) -> str:
    return str(value or "-").replace("|", "\\|").replace("\n", " ")


def _render_ghs_sources(classifications: Sequence[Dict[str, Any]]) -> str:
    rendered: List[str] = []
    for item in classifications:
        source_name = str(item.get("source_name") or "unnamed source")
        source_url = str(item.get("source_url") or "")
        reference = item.get("reference_number")
        label = f"[{source_name}]({source_url})" if source_url else source_name
        if reference:
            label += f" (ref {reference})"
        statements = list(item.get("hazard_statements") or [])
        detail = "; ".join(statements[:4]) or "classification present"
        if len(statements) > 4:
            detail += f"; +{len(statements) - 4} more"
        license_url = str(item.get("license_url") or "")
        if license_url:
            detail += f"; [license]({license_url})"
        elif item.get("license_note"):
            detail += f"; license: {item.get('license_note')}"
        rendered.append(f"{label}: {detail}")
    return "<br/>".join(rendered) if rendered else "-"


def render_post_route_audit_markdown(
    payload: Dict[str, Any],
    *,
    heading_level: int = 1,
) -> str:
    level = min(max(int(heading_level), 1), 5)
    heading = "#" * level
    section_heading = "#" * (level + 1)
    incomplete_text = ", ".join(payload.get("incompleteness_reasons") or [])
    lines = [
        f"{heading} Rachel Post-Finalize Audit",
        "",
        f"- audit id: `{payload.get('audit_id', '')}`",
        f"- route digest: `{payload.get('route_digest', '')}`",
        f"- generated at: {payload.get('generated_at', '')}",
        f"- status: **{payload.get('status', 'unknown')}**",
        f"- incompleteness reasons: {_md_cell(incomplete_text)}",
        "",
    ]
    buyability = payload.get("buyability") or {}
    if buyability:
        lines.extend([
            f"{section_heading} Buyability",
            "",
            "PubChem B1/B2 evidence is identity/catalog evidence only; it is not live stock, price, lead time, or procurement confirmation.",
            "",
            "| node | SMILES | grade | CID | vendor catalog | example vendor records |",
            "|---|---|---:|---:|---|---|",
        ])
        for item in buyability.get("terminals", []) or []:
            catalog = item.get("vendor_catalog") or {}
            examples = []
            for row in list(catalog.get("rows") or [])[:3]:
                label = f"{row.get('source_name', '')} {row.get('registry_id', '')}".strip()
                url = row.get("source_record_url") or row.get("source_url")
                examples.append(f"[{label}]({url})" if url else label)
            lines.append(
                "| {node} | `{smiles}` | {grade} | {cid} | {status} | {vendors} |".format(
                    node=_md_cell(item.get("node_id")),
                    smiles=_md_cell(item.get("smiles")),
                    grade=_md_cell(item.get("evidence_grade")),
                    cid=_md_cell((item.get("pubchem") or {}).get("best_cid")),
                    status=_md_cell(catalog.get("status")),
                    vendors="<br/>".join(examples) if examples else "-",
                )
            )
        lines.append("")

    safety = payload.get("safety") or {}
    if safety:
        lines.extend([
            f"{section_heading} Substance Safety Evidence",
            "",
            "Each GHS classification below remains attached to its PubChem depositor/authority source and license; conflicting sources are not collapsed into a single verdict.",
            "",
            "| material | roles | CID | grade | checked-source state | local structural alerts | source-attributed GHS evidence |",
            "|---|---|---:|---:|---|---|---|",
        ])
        for item in safety.get("substances", []) or []:
            ghs = item.get("ghs") or {}
            structural_alerts = "; ".join(
                str(alert.get("alert_id") or "")
                for alert in item.get("structural_screening", []) or []
                if alert.get("alert_id")
            )
            lines.append(
                "| `{smiles}` | {roles} | {cid} | {grade} | {state} | {structural} | {sources} |".format(
                    smiles=_md_cell(item.get("smiles")),
                    roles=_md_cell(", ".join(item.get("roles") or [])),
                    cid=_md_cell(item.get("cid")),
                    grade=_md_cell(item.get("evidence_grade")),
                    state=_md_cell(item.get("status")),
                    structural=_md_cell(structural_alerts),
                    sources=_md_cell(_render_ghs_sources(ghs.get("classifications") or [])),
                )
            )
        lines.extend([
            "",
            f"{section_heading} Step Safety Screening",
            "",
            "The explicit assessment is LLM-authored planning judgment; its semantic correctness is not independently validated.",
            "",
            "| step | reaction | explicit assessment | local alerts | condition gaps | process status |",
            "|---|---|---|---|---|---|",
        ])
        for item in safety.get("steps", []) or []:
            alerts = "; ".join(
                dangerous_combo_text(
                    str((alert.get("groups") or [""])[0]),
                    str(alert.get("rationale", "") or ""),
                )
                for alert in item.get("alerts", [])
            )
            lines.append(
                "| {step} | {reaction} | {declared} | {alerts} | {gaps} | {status} |".format(
                    step=_md_cell(f"Step {item.get('forward_step', '?')} ({item.get('step_id', '')})"),
                    reaction=_md_cell(item.get("reaction_name")),
                    declared=_md_cell(item.get("declared_risk_assessment") or "not assessed"),
                    alerts=_md_cell(alerts),
                    gaps=_md_cell(format_process_condition_gaps(item.get("process_condition_coverage") or {})),
                    status=_md_cell(item.get("process_assessment_status")),
                )
            )
        lines.append("")

    lines.extend([f"{section_heading} Limitations", ""])
    limitations = list(payload.get("limitations") or [])
    limitations.extend(buyability.get("limitations") or [])
    limitations.extend(safety.get("limitations") or [])
    for limitation in _unique_strings(limitations):
        lines.append(f"- {limitation}")
    return "\n".join(lines) + "\n"


def render_post_route_audit_html(payload: Dict[str, Any]) -> str:
    """Render the full current audit inside the self-contained HTML report."""
    audit_markdown = render_post_route_audit_markdown(payload, heading_level=3)
    return (
        '<div class="section post-route-audit">'
        '<h2>Post-Finalize Safety and Buyability Audit</h2>'
        '<p>This independent evidence snapshot does not change route chemistry, '
        'validation gates, terminal decisions, or ranking.</p>'
        f'<pre class="tree-text">{html.escape(audit_markdown)}</pre>'
        '</div>'
    )


def write_post_route_audit_files(payload: Dict[str, Any], output_dir: Path) -> List[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "post_route_audit.json"
    md_path = output_dir / "post_route_audit.md"
    json_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    md_path.write_text(render_post_route_audit_markdown(payload), encoding="utf-8")
    return [str(json_path), str(md_path)]


def render_post_route_audit_summary(payload: Dict[str, Any]) -> str:
    if not payload:
        return ""
    lines = ["", "Post-finalize audit", "-" * 40]
    lines.append(f"Status: {payload.get('status', 'unknown')}")
    incomplete_reasons = payload.get("incompleteness_reasons") or []
    if incomplete_reasons:
        lines.append("Incomplete because: " + ", ".join(incomplete_reasons))
    buyability = payload.get("buyability") or {}
    if buyability:
        grades = (buyability.get("summary") or {}).get("by_evidence_grade") or {}
        lines.append(
            "Buyability evidence: "
            + ", ".join(f"{grade}={grades.get(grade, 0)}" for grade in ("B0", "B1", "B2"))
        )
    safety = payload.get("safety") or {}
    if safety:
        summary = safety.get("summary") or {}
        lines.append(
            "Safety evidence: "
            f"S1 substances={summary.get('substances_with_source_attributed_ghs', 0)}, "
            f"steps missing explicit assessment={summary.get('steps_missing_explicit_risk_assessment', 0)}"
        )
        ready = summary.get("steps_with_process_conditions_for_preliminary_review", 0)
        incomplete = summary.get("steps_with_incomplete_process_conditions", 0)
        lines.append(
            "Process condition coverage: "
            f"preliminary-review-ready={ready}, incomplete={incomplete}; "
            "field coverage is not a safe/unsafe verdict."
        )
    return "\n".join(lines)
