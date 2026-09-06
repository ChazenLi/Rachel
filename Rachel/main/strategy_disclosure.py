"""Site-first and reaction-first action disclosure for Rachel sessions.

This module normalizes bond, smart-capping, FGI, and custom precursor proposals
into stable internal candidate IDs. The LLM-facing path exposes those proposal
units as actions: site_reaction_map -> explore_site -> try_action.
explore_reaction remains an auxiliary view.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional


SMART_CAPPING_PUBLIC_LABEL = "LLM-completed structural disconnection"


@dataclass
class CandidateUnit:
    """Normalized candidate wrapper across bond, smart-cap, FGI, and terminal options."""

    candidate_id: str
    source: str
    source_label: str
    reaction_id: str = ""
    reaction_type: str = ""
    template_id: str = ""
    template_name: str = ""
    precursors_preview: List[str] = field(default_factory=list)
    reagents: List[str] = field(default_factory=list)
    fragments: List[str] = field(default_factory=list)
    bond_idx: Optional[int] = None
    alt_idx: Optional[int] = None
    fgi_idx: Optional[int] = None
    smart_idx: Optional[int] = None
    actual_bond_idx: Optional[int] = None
    atoms: List[int] = field(default_factory=list)
    role_pair: List[str] = field(default_factory=list)
    bond_type: str = ""
    in_ring: bool = False
    heuristic_score: float = 0.0
    confidence: Optional[float] = None
    description: str = ""
    candidate_label: str = ""
    why_existing_candidates_rejected: str = ""
    rationale_summary: str = ""
    intended_deltas: List[str] = field(default_factory=list)
    expected_ring_change: str = ""
    changed_bonds: List[Dict[str, Any]] = field(default_factory=list)
    preserved_anchors: List[Any] = field(default_factory=list)
    mechanistic_evidence: List[str] = field(default_factory=list)
    family_evidence: Dict[str, Any] = field(default_factory=dict)
    duplicate_of: str = ""
    risk_tags: List[str] = field(default_factory=list)
    experience_card_hints: List[str] = field(default_factory=list)
    strategy_id: str = ""
    route_sketch_id: str = ""
    route_strategy_brief: Dict[str, Any] = field(default_factory=dict)
    rescue_id: str = ""
    rescue_step_idx: Optional[int] = None
    continuation_precursor: str = ""
    precursor_normalization: Optional[Dict[str, Any]] = None
    execution: Dict[str, Any] = field(default_factory=dict)
    knowledge_refs: List[Dict[str, str]] = field(default_factory=list)

    def to_dict(self, *, compact: bool = False) -> Dict[str, Any]:
        data = asdict(self)
        if compact:
            keep = {
                "candidate_id",
                "source",
                "source_label",
                "reaction_id",
                "reaction_type",
                "template_id",
                "template_name",
                "precursors_preview",
                "reagents",
                "bond_idx",
                "alt_idx",
                "fgi_idx",
                "smart_idx",
                "actual_bond_idx",
                "atoms",
                "role_pair",
                "heuristic_score",
                "confidence",
                "candidate_label",
                "intended_deltas",
                "expected_ring_change",
                "changed_bonds",
                "preserved_anchors",
                "mechanistic_evidence",
                "family_evidence",
                "duplicate_of",
                "risk_tags",
                "route_sketch_id",
                "route_strategy_brief",
                "rescue_id",
                "rescue_step_idx",
                "continuation_precursor",
                "precursor_normalization",
                "execution",
                "knowledge_refs",
            }
            data = {k: v for k, v in data.items() if k in keep}
        return {k: v for k, v in data.items() if v not in (None, "", [], {})}


def _contains_attachment_placeholder(values: Iterable[Any]) -> bool:
    """Return whether a precursor preview contains an intentional dummy atom."""
    return any("*" in str(value) for value in (values or []))


def _requires_llm_completion(candidate: CandidateUnit) -> bool:
    return "requires_llm_completion" in candidate.risk_tags


def _completion_obligation() -> Dict[str, Any]:
    return {
        "task": "construct a complete one-step reaction proposal",
        "check": [
            "chemically correct closed-shell precursors or justified in-situ species",
            "missing reagents or atom-source components",
            "mechanism, chemoselectivity, and functional-group compatibility",
        ],
        "next_step": "propose_action",
    }


def _candidate_payload_to_action(data: Dict[str, Any]) -> Dict[str, Any]:
    """Project an internal candidate payload into the LLM-facing action schema."""
    replacements = {
        "candidate_id": "action_id",
        "candidate_label": "action_label",
        "why_existing_candidates_rejected": "why_existing_actions_rejected",
    }
    action: Dict[str, Any] = {}
    for key, value in data.items():
        action[replacements.get(key, key)] = value
    return action


def _action_dict(candidate: CandidateUnit, *, compact: bool = False) -> Dict[str, Any]:
    return _candidate_payload_to_action(candidate.to_dict(compact=compact))


def _site_aligned_action_id(candidate: CandidateUnit) -> str:
    """Return the public action handle aligned with the real site id.

    Earlier bond candidates used their position in the disconnectable-bond list
    as the first numeric component. The first-layer site menu uses
    actual_bond_idx, so expose the same index in second-layer action ids while
    keeping legacy ids accepted by lookup.
    """
    if candidate.source == "bond" and candidate.alt_idx is not None:
        bond_ref = candidate.actual_bond_idx
        if bond_ref is None:
            bond_ref = candidate.bond_idx
        if bond_ref is not None:
            return f"site:bond:{bond_ref}:alt:{candidate.alt_idx}"
    if candidate.source == "smart_capping" and candidate.smart_idx is not None:
        bond_ref = candidate.actual_bond_idx
        if bond_ref is None:
            bond_ref = candidate.bond_idx
        if bond_ref is not None:
            return f"site:bond:{bond_ref}:smart:{candidate.smart_idx}"
    return candidate.candidate_id


def _legacy_action_ids(candidate: CandidateUnit) -> List[str]:
    aliases: List[str] = []
    if candidate.source == "bond" and candidate.bond_idx is not None and candidate.alt_idx is not None:
        aliases.append(f"bond:{candidate.bond_idx}:alt:{candidate.alt_idx}")
    if candidate.source == "smart_capping" and candidate.bond_idx is not None and candidate.smart_idx is not None:
        aliases.append(f"smart:{candidate.bond_idx}:{candidate.smart_idx}")
    return _unique_in_order(alias for alias in aliases if alias and alias != candidate.candidate_id)


def _reaction_type(template: str) -> str:
    if not template:
        return ""
    return _clean_reaction_name(template.split("(")[0].strip())


def _clean_reaction_name(name: str) -> str:
    """Normalize reaction names for LLM-facing disclosure."""
    if not name:
        return ""
    cleaned = str(name)
    if "\u00e2\u20ac" in cleaned:
        try:
            cleaned = cleaned.encode("cp1252").decode("utf-8")
        except UnicodeError:
            pass
    cleaned = (
        cleaned.replace("\u00e2\u20ac\u201c", "-")
        .replace("\u00e2\u20ac\u201d", "-")
        .replace("\u2010", "-")
        .replace("\u2011", "-")
        .replace("\u2012", "-")
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2015", "-")
        .replace("\ufffdC", "-")
        .replace("\ufffd", "-")
    )
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned


def _list_field(value: Any) -> List[Any]:
    """Normalize optional list-like public payload fields without splitting strings."""
    if value in (None, "", [], {}):
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def _reaction_id(reaction_name: str) -> str:
    base = _clean_reaction_name(reaction_name).lower()
    slug = re.sub(r"[^a-z0-9]+", "_", base).strip("_")
    return f"reaction:{slug or 'unknown'}"


def _candidate_reaction_name(candidate: CandidateUnit) -> str:
    return _clean_reaction_name(
        candidate.reaction_type or candidate.template_name or candidate.source
    )


def _candidate_reaction_id(candidate: CandidateUnit) -> str:
    return candidate.reaction_id or _reaction_id(_candidate_reaction_name(candidate))


def _risk_tags(
    *,
    source: str,
    bond: Optional[Dict[str, Any]] = None,
    payload: Optional[Dict[str, Any]] = None,
) -> List[str]:
    bond = bond or {}
    payload = payload or {}
    tags: List[str] = []
    if bond.get("in_ring"):
        tags.append("ring_bond")
    if payload.get("incompatible_with"):
        tags.append("functional_group_compatibility")
    if not (payload.get("precursors") or payload.get("fragments")) and source != "terminal":
        tags.append("no_precursor_preview")
    preview = payload.get("precursors") or payload.get("fragments") or []
    if source in {"bond", "fgi", "smart_capping"} and _contains_attachment_placeholder(preview):
        tags.append("intentional_attachment_placeholder")
        tags.append("requires_llm_completion")
    if source == "smart_capping":
        tags.append("heuristic_cap")
        if "requires_llm_completion" not in tags:
            tags.append("requires_llm_completion")
    if source == "custom_precursors":
        tags.append("llm_custom_precursor")
    if source == "terminal":
        tags.append("requires_buyability_or_simplicity_check")
    return tags


def _build_bond_candidate(
    bond: Dict[str, Any],
    bond_idx: int,
    alt: Dict[str, Any],
    alt_idx: int,
) -> CandidateUnit:
    reaction = _reaction_type(alt.get("template", ""))
    precursors = list(alt.get("precursors", []) or [])
    atoms = list(bond.get("atoms", []) or [])
    risk_tags = _risk_tags(source="bond", bond=bond, payload=alt)
    execution: Dict[str, Any] = {
        "mode": "retro_template",
        "template_id": alt.get("template_id", ""),
        "bond": atoms,
        "expected_precursors": precursors,
    }
    if "intentional_attachment_placeholder" in risk_tags:
        execution = {
            "mode": "llm_completion_required",
            "template_id": alt.get("template_id", ""),
            "bond": atoms,
            "structural_hint_precursors": precursors,
            "placeholder_semantics": (
                "The dummy atom marks an intentional retrosynthetic attachment site; "
                "the preview is not a real precursor set or a template-failure verdict."
            ),
            "completion_obligation": _completion_obligation(),
        }
    return CandidateUnit(
        candidate_id=f"bond:{bond_idx}:alt:{alt_idx}",
        source="bond",
        source_label=f"bond[{bond_idx}].alt[{alt_idx}]",
        reaction_type=reaction,
        template_id=alt.get("template_id", ""),
        template_name=alt.get("template", ""),
        precursors_preview=precursors,
        bond_idx=bond_idx,
        alt_idx=alt_idx,
        actual_bond_idx=bond.get("actual_bond_idx"),
        atoms=atoms,
        role_pair=list(bond.get("role_pair", []) or []),
        bond_type=bond.get("bond_type", ""),
        in_ring=bool(bond.get("in_ring", False)),
        heuristic_score=float(bond.get("heuristic_score", 0) or 0),
        risk_tags=risk_tags,
        execution=execution,
        knowledge_refs=[alt["knowledge_ref"]]
        if alt.get("knowledge_ref")
        else [],
    )


def _build_smart_candidate(
    bond: Dict[str, Any],
    bond_idx: int,
    cap: Dict[str, Any],
    smart_idx: int,
) -> CandidateUnit:
    fragments = list(cap.get("fragments", []) or [])
    heuristic_reaction_hint = _clean_reaction_name(cap.get("reaction_type", ""))
    return CandidateUnit(
        candidate_id=f"smart:{bond_idx}:{smart_idx}",
        source="smart_capping",
        source_label=f"bond[{bond_idx}].smart_capping[{smart_idx}]",
        reaction_type=SMART_CAPPING_PUBLIC_LABEL,
        precursors_preview=fragments,
        fragments=fragments,
        bond_idx=bond_idx,
        smart_idx=smart_idx,
        actual_bond_idx=bond.get("actual_bond_idx"),
        atoms=list(bond.get("atoms", []) or []),
        role_pair=list(bond.get("role_pair", []) or []),
        bond_type=bond.get("bond_type", ""),
        in_ring=bool(bond.get("in_ring", False)),
        heuristic_score=float(bond.get("heuristic_score", 0) or 0),
        confidence=cap.get("confidence"),
        description=cap.get("description", ""),
        risk_tags=_risk_tags(source="smart_capping", bond=bond, payload=cap),
        execution={
            "mode": "llm_completion_required",
            "heuristic_reaction_hint": heuristic_reaction_hint,
            "structural_hint_precursors": fragments,
            "completion_obligation": _completion_obligation(),
        },
        knowledge_refs=[cap["knowledge_ref"]]
        if cap.get("knowledge_ref")
        else [],
    )


def _build_fgi_candidate(fgi: Dict[str, Any], fgi_idx: int) -> CandidateUnit:
    return CandidateUnit(
        candidate_id=f"fgi:{fgi_idx}",
        source="fgi",
        source_label=f"fgi[{fgi_idx}]",
        reaction_type=_reaction_type(fgi.get("template", "")) or "fgi",
        template_id=fgi.get("template_id", ""),
        template_name=fgi.get("template", ""),
        precursors_preview=list(fgi.get("precursors", []) or []),
        fgi_idx=fgi_idx,
        risk_tags=_risk_tags(source="fgi", payload=fgi),
        execution={
            "mode": "fgi_template",
            "template_id": fgi.get("template_id", ""),
            "expected_precursors": list(fgi.get("precursors", []) or []),
        },
        knowledge_refs=[fgi["knowledge_ref"]]
        if fgi.get("knowledge_ref")
        else [],
    )


def _terminal_candidate() -> CandidateUnit:
    return CandidateUnit(
        candidate_id="terminal:accept",
        source="terminal",
        source_label="accept_terminal",
        reaction_type="Terminal Acceptance",
        risk_tags=_risk_tags(source="terminal"),
        execution={"mode": "terminal_acceptance"},
    )


def _build_custom_candidate(custom: Dict[str, Any], custom_idx: int) -> CandidateUnit:
    strategy_id = str(custom.get("strategy_id", "") or "").strip()
    slug_seed = strategy_id or "llm_custom"
    custom_slug = re.sub(r"[^A-Za-z0-9]+", "_", str(slug_seed)).strip("_").lower() or "llm_custom"
    precursors = list(custom.get("precursors", []) or custom.get("precursors_preview", []) or [])
    reaction_name = _clean_reaction_name(
        custom.get("reaction_name") or custom.get("reaction_type") or "llm_proposed"
    )
    explicit_reaction_id = custom.get("reaction_id") or _reaction_id(reaction_name)
    risk_tags = list(custom.get("risk_tags", []) or [])
    for tag in _risk_tags(source="custom_precursors", payload={"precursors": precursors}):
        if tag not in risk_tags:
            risk_tags.append(tag)
    return CandidateUnit(
        candidate_id=custom.get("candidate_id", f"custom:{custom_slug}:{custom_idx}"),
        source="custom_precursors",
        source_label=custom.get("source_label", f"custom[{custom_slug}:{custom_idx}]"),
        reaction_id=explicit_reaction_id,
        reaction_type=reaction_name,
        template_name=reaction_name,
        precursors_preview=precursors,
        reagents=list(custom.get("reagents", []) or []),
        fragments=precursors,
        description=custom.get("rationale_summary", ""),
        candidate_label=custom.get("candidate_label", ""),
        why_existing_candidates_rejected=custom.get("why_existing_candidates_rejected", ""),
        rationale_summary=custom.get("rationale_summary", ""),
        intended_deltas=[str(item) for item in _list_field(custom.get("intended_deltas"))],
        expected_ring_change=custom.get("expected_ring_change", ""),
        changed_bonds=[
            dict(item) if isinstance(item, dict) else {"description": item}
            for item in _list_field(custom.get("changed_bonds"))
        ],
        preserved_anchors=_list_field(custom.get("preserved_anchors")),
        mechanistic_evidence=[str(item) for item in _list_field(custom.get("mechanistic_evidence"))],
        family_evidence=dict(custom.get("family_evidence", {}) or {}),
        in_ring=custom.get("in_ring") if "in_ring" in custom else None,
        duplicate_of=custom.get("duplicate_of", ""),
        risk_tags=risk_tags,
        experience_card_hints=list(custom.get("experience_card_hints", []) or []),
        strategy_id=strategy_id,
        route_sketch_id=custom.get("route_sketch_id", ""),
        route_strategy_brief=dict(custom.get("route_strategy_brief", {}) or {}),
        rescue_id=custom.get("rescue_id", ""),
        rescue_step_idx=custom.get("rescue_step_idx"),
        continuation_precursor=custom.get("continuation_precursor", ""),
        precursor_normalization=custom.get("precursor_normalization"),
        execution={
            "mode": "custom_precursors",
            "expected_precursors": precursors,
            "reagents": list(custom.get("reagents", []) or []),
        },
    )


def build_candidate_units(decision_context: Dict[str, Any]) -> List[CandidateUnit]:
    """Normalize all currently generated proposal candidates."""
    units: List[CandidateUnit] = []
    for bond_idx, bond in enumerate(decision_context.get("disconnectable_bonds", []) or []):
        for alt_idx, alt in enumerate(bond.get("alternatives", []) or []):
            units.append(_build_bond_candidate(bond, bond_idx, alt, alt_idx))
        for smart_idx, cap in enumerate(bond.get("smart_capping", []) or []):
            units.append(_build_smart_candidate(bond, bond_idx, cap, smart_idx))

    for fgi_idx, fgi in enumerate(decision_context.get("fgi_options", []) or []):
        units.append(_build_fgi_candidate(fgi, fgi_idx))

    for custom_idx, custom in enumerate(decision_context.get("custom_candidates", []) or []):
        units.append(_build_custom_candidate(custom, custom_idx))

    units.append(_terminal_candidate())
    return units


def terminal_acceptance_action() -> Dict[str, Any]:
    """Expose terminal acceptance as a route action, not as a reaction family."""
    terminal = _terminal_candidate().to_dict(compact=True)
    terminal["next_step"] = "accept_terminal(reason=...)"
    terminal["recommended_check"] = "Confirm buyability, simplicity, or route stopping criteria before accepting."
    return terminal


def _merge_tags(candidates: List[CandidateUnit]) -> List[str]:
    seen = set()
    tags: List[str] = []
    for candidate in candidates:
        for tag in candidate.risk_tags:
            if tag not in seen:
                seen.add(tag)
                tags.append(tag)
    return tags


def _candidate_site_hint(candidate: CandidateUnit) -> str:
    if candidate.source == "fgi":
        label = candidate.reaction_type or candidate.template_name or "functional-group option"
        return f"FGI: {_clean_reaction_name(label)}"

    if candidate.source == "custom_precursors":
        label = candidate.reaction_type or candidate.candidate_label or "custom precursors"
        return f"custom: {_clean_reaction_name(label)}"

    if candidate.role_pair:
        bond_label = "-".join(str(role) for role in candidate.role_pair if role)
        if bond_label:
            suffix = "ring bond" if candidate.in_ring else "bond"
            return f"{bond_label} {suffix}"

    if candidate.atoms:
        atom_label = "-".join(str(atom) for atom in candidate.atoms)
        suffix = "ring bond" if candidate.in_ring else "bond"
        return f"atoms {atom_label} {suffix}"

    return candidate.source or "candidate"


def _reaction_site_hint(candidates: List[CandidateUnit]) -> str:
    hints = _unique_in_order(_candidate_site_hint(candidate) for candidate in candidates)
    if not hints:
        return ""
    if len(hints) == 1:
        return hints[0]
    if len(hints) <= 3:
        return "; ".join(hints)
    return f"{len(hints)} sites: " + "; ".join(hints[:3])


def _risk_hint(candidates: List[CandidateUnit]) -> str:
    tags = _merge_tags(candidates)
    return ", ".join(tags) if tags else "none"


def _candidate_site_id(candidate: CandidateUnit) -> str:
    if candidate.source in {"bond", "smart_capping"}:
        if candidate.actual_bond_idx is not None:
            return f"bond:{candidate.actual_bond_idx}"
        if candidate.bond_idx is not None:
            return f"bond:{candidate.bond_idx}"
    if candidate.source == "fgi":
        return f"fgi:{candidate.fgi_idx if candidate.fgi_idx is not None else 'unknown'}"
    if candidate.source == "custom_precursors":
        return f"custom:{candidate.candidate_id}"
    return f"{candidate.source or 'candidate'}:{candidate.candidate_id}"


def _candidate_site_type(candidate: CandidateUnit) -> str:
    if candidate.source in {"bond", "smart_capping"}:
        return "bond"
    if candidate.source == "fgi":
        return "fgi"
    if candidate.source == "custom_precursors":
        return "custom"
    return candidate.source or "candidate"


def _unique_in_order(values: Iterable[str]) -> List[str]:
    seen = set()
    result: List[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def build_reaction_index(
    decision_context: Dict[str, Any],
    *,
    representative_limit: int = 2,
) -> List[Dict[str, Any]]:
    """Build an internal diagnostic index grouped by concrete reaction name."""
    units = build_candidate_units(decision_context)
    by_reaction: Dict[str, List[CandidateUnit]] = {}
    reaction_names: Dict[str, str] = {}

    for unit in units:
        if unit.source == "terminal":
            continue
        reaction_name = _clean_reaction_name(unit.reaction_type or unit.template_name or unit.source)
        rid = unit.reaction_id or _reaction_id(reaction_name)
        by_reaction.setdefault(rid, []).append(unit)
        reaction_names.setdefault(rid, reaction_name)

    families: List[Dict[str, Any]] = []
    for rid, candidates in by_reaction.items():
        representatives = candidates[: max(0, representative_limit)]
        sources = _unique_in_order(candidate.source for candidate in candidates)
        families.append(
            {
                "reaction_id": rid,
                "reaction_name": reaction_names[rid],
                "sources": sources,
                "candidate_count": len(candidates),
                "candidate_ids": [candidate.candidate_id for candidate in candidates],
                "representative_candidates": [
                    candidate.to_dict(compact=True) for candidate in representatives
                ],
                "risk_tags": _merge_tags(candidates),
                "site_hint": _reaction_site_hint(candidates),
                "risk_hint": _risk_hint(candidates),
                "next_step": f"explore_reaction('{rid}')",
            }
        )
    return families


def build_reaction_menu(decision_context: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build the first-layer LLM menu: all reaction names, no action detail."""
    menu: List[Dict[str, Any]] = []
    for family in build_reaction_index(decision_context):
        menu.append(
            {
                "reaction_id": family["reaction_id"],
                "reaction_name": family["reaction_name"],
                "action_count": family["candidate_count"],
                "source_summary": family["sources"],
                "site_hint": family["site_hint"],
                "risk_hint": family["risk_hint"],
            }
        )
    return menu


def _site_reaction_row(reaction_id: str, reaction_name: str, candidates: List[CandidateUnit]) -> Dict[str, Any]:
    return {
        "reaction_id": reaction_id,
        "reaction_name": reaction_name,
        "action_count": len(candidates),
        "source_summary": _unique_in_order(candidate.source for candidate in candidates),
        "risk_hint": _risk_hint(candidates),
    }


def _candidate_site_detail(candidate: CandidateUnit) -> Dict[str, Any]:
    """Return the LLM-facing second-layer action view for one site."""
    reaction_name = _candidate_reaction_name(candidate)
    data: Dict[str, Any] = {
        "action_id": _site_aligned_action_id(candidate),
        "site_id": _candidate_site_id(candidate),
        "reaction": {
            "id": _candidate_reaction_id(candidate),
            "name": reaction_name,
        },
        "source": candidate.source,
        "source_label": candidate.source_label,
        "template_id": candidate.template_id,
        "template_name": candidate.template_name,
        "precursors_preview": list(candidate.precursors_preview),
        "reagents": list(candidate.reagents),
        "bond_idx": candidate.bond_idx,
        "alt_idx": candidate.alt_idx,
        "fgi_idx": candidate.fgi_idx,
        "smart_idx": candidate.smart_idx,
        "actual_bond_idx": candidate.actual_bond_idx,
        "atoms": list(candidate.atoms),
        "role_pair": list(candidate.role_pair),
        "bond_type": candidate.bond_type,
        "in_ring": candidate.in_ring,
        "heuristic_score": candidate.heuristic_score,
        "confidence": candidate.confidence,
        "action_label": candidate.candidate_label,
        "why_existing_actions_rejected": candidate.why_existing_candidates_rejected,
        "rationale_summary": candidate.rationale_summary,
        "intended_deltas": list(candidate.intended_deltas),
        "expected_ring_change": candidate.expected_ring_change,
        "changed_bonds": list(candidate.changed_bonds),
        "preserved_anchors": list(candidate.preserved_anchors),
        "mechanistic_evidence": list(candidate.mechanistic_evidence),
        "family_evidence": dict(candidate.family_evidence),
        "duplicate_of": candidate.duplicate_of,
        "risk_tags": list(candidate.risk_tags),
        "execution": dict(candidate.execution),
        "knowledge_refs": list(candidate.knowledge_refs),
        "next_step": (
            "propose_action"
            if _requires_llm_completion(candidate)
            else "try_action(action_id)"
        ),
    }
    return {k: v for k, v in data.items() if v not in (None, "", [], {})}


def build_site_reaction_map(decision_context: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Group compact reactions by the real bond/FGI/custom site they act on."""
    units = [unit for unit in build_candidate_units(decision_context) if unit.source != "terminal"]
    by_site: Dict[str, List[CandidateUnit]] = {}
    site_order: List[str] = []

    for unit in units:
        site_id = _candidate_site_id(unit)
        if site_id not in by_site:
            by_site[site_id] = []
            site_order.append(site_id)
        by_site[site_id].append(unit)

    site_map: List[Dict[str, Any]] = []
    for site_id in site_order:
        site_candidates = by_site[site_id]
        by_reaction: Dict[str, List[CandidateUnit]] = {}
        reaction_names: Dict[str, str] = {}
        reaction_order: List[str] = []
        for candidate in site_candidates:
            reaction_name = _clean_reaction_name(
                candidate.reaction_type or candidate.template_name or candidate.source
            )
            reaction_id = candidate.reaction_id or _reaction_id(reaction_name)
            if reaction_id not in by_reaction:
                by_reaction[reaction_id] = []
                reaction_names[reaction_id] = reaction_name
                reaction_order.append(reaction_id)
            by_reaction[reaction_id].append(candidate)

        reactions = [
            _site_reaction_row(reaction_id, reaction_names[reaction_id], by_reaction[reaction_id])
            for reaction_id in reaction_order
        ]
        site_payload = {
            "site_id": site_id,
            "site_type": _candidate_site_type(site_candidates[0]),
            "site_hint": _reaction_site_hint(site_candidates),
            "action_count": len(site_candidates),
            "reaction_count": len(reactions),
            "source_summary": _unique_in_order(candidate.source for candidate in site_candidates),
            "risk_hint": _risk_hint(site_candidates),
            "reactions": reactions,
            "next_step": f"explore_site('{site_id}')",
        }
        if len(reactions) > 1:
            site_payload["competition_hint"] = f"same-site competition: {len(reactions)} reactions"
        site_map.append(site_payload)
    return site_map


def build_reaction_opportunity_brief(
    decision_context: Dict[str, Any],
    *,
    high_competition_limit: int = 3,
) -> Dict[str, Any]:
    """Summarize reaction opportunities without embedding the full first layer."""
    site_map = build_site_reaction_map(decision_context)
    reaction_names = _unique_in_order(
        reaction.get("reaction_name", "")
        for site in site_map
        for reaction in site.get("reactions", []) or []
    )
    competing_sites = [
        site for site in site_map if int(site.get("reaction_count", 0) or 0) > 1
    ]
    ranked_competing_sites = sorted(
        competing_sites,
        key=lambda site: (
            int(site.get("reaction_count", 0) or 0),
            int(site.get("action_count", 0) or 0),
        ),
        reverse=True,
    )
    high_competition_sites: List[Dict[str, Any]] = []
    for site in ranked_competing_sites[: max(0, high_competition_limit)]:
        high_competition_sites.append(
            {
                "site_id": site.get("site_id", ""),
                "site_hint": site.get("site_hint", ""),
                "reaction_count": site.get("reaction_count", 0),
                "action_count": site.get("action_count", 0),
                "reaction_names": [
                    reaction.get("reaction_name", "")
                    for reaction in site.get("reactions", []) or []
                    if reaction.get("reaction_name", "")
                ],
                "risk_hint": site.get("risk_hint", "none"),
                "next_step": site.get("next_step", ""),
            }
        )

    return {
        "site_count": len(site_map),
        "total_reaction_count": sum(
            len(site.get("reactions", []) or []) for site in site_map
        ),
        "action_count": sum(
            int(site.get("action_count", 0) or 0) for site in site_map
        ),
        "competing_site_count": len(competing_sites),
        "high_competition_sites": high_competition_sites,
        "reaction_names": reaction_names,
        "first_layer_command": "reaction_sites()",
    }


def explore_site(decision_context: Dict[str, Any], site_id: str) -> Dict[str, Any]:
    """Return all candidates competing at one real bond/FGI/custom site."""
    units = [
        unit
        for unit in build_candidate_units(decision_context)
        if unit.source != "terminal"
    ]
    matching = [unit for unit in units if _candidate_site_id(unit) == site_id]
    if not matching:
        available = [
            {
                "site_id": site.get("site_id", ""),
                "site_hint": site.get("site_hint", ""),
                "action_count": site.get("action_count", 0),
                "reaction_count": site.get("reaction_count", 0),
            }
            for site in build_site_reaction_map(decision_context)
        ]
        return {"error": f"unknown site_id: {site_id}", "available": available}

    by_reaction: Dict[str, List[CandidateUnit]] = {}
    reaction_names: Dict[str, str] = {}
    reaction_order: List[str] = []
    for candidate in matching:
        rid = _candidate_reaction_id(candidate)
        if rid not in by_reaction:
            by_reaction[rid] = []
            reaction_names[rid] = _candidate_reaction_name(candidate)
            reaction_order.append(rid)
        by_reaction[rid].append(candidate)

    first = matching[0]
    executable_actions = [candidate for candidate in matching if not _requires_llm_completion(candidate)]
    placeholder_completion_required = any(
        "intentional_attachment_placeholder" in candidate.risk_tags
        for candidate in matching
    )
    result: Dict[str, Any] = {
        "site_id": site_id,
        "site_type": _candidate_site_type(first),
        "site_hint": _reaction_site_hint(matching),
        "action_count": len(matching),
        "reaction_count": len(reaction_order),
        "source_summary": _unique_in_order(candidate.source for candidate in matching),
        "risk_tags": _merge_tags(matching),
        "risk_hint": _risk_hint(matching),
        "reactions": [
            _site_reaction_row(rid, reaction_names[rid], by_reaction[rid])
            for rid in reaction_order
        ],
        "actions": [_candidate_site_detail(candidate) for candidate in matching],
        "next_step": (
            "try_action(action_id)"
            if executable_actions
            else "propose_action"
        ),
    }
    if placeholder_completion_required:
        result["open_action_requires_completion"] = True
    if len(reaction_order) > 1:
        result["competition_hint"] = f"same-site competition: {len(reaction_order)} reactions"
    for key in ("actual_bond_idx", "atoms", "role_pair", "bond_type", "in_ring", "heuristic_score"):
        value = getattr(first, key)
        if value not in (None, "", [], {}):
            result[key] = value
    return result


def _same_site_context(
    units: List[CandidateUnit],
    matching: List[CandidateUnit],
    current_reaction_id: str,
) -> Dict[str, Any]:
    matching_site_ids = _unique_in_order(_candidate_site_id(candidate) for candidate in matching)
    if len(matching_site_ids) != 1:
        return {}

    site_id = matching_site_ids[0]
    site_candidates = [
        candidate
        for candidate in units
        if candidate.source != "terminal" and _candidate_site_id(candidate) == site_id
    ]
    by_reaction: Dict[str, str] = {}
    reaction_order: List[str] = []
    for candidate in site_candidates:
        rid = _candidate_reaction_id(candidate)
        if rid not in by_reaction:
            by_reaction[rid] = _candidate_reaction_name(candidate)
            reaction_order.append(rid)

    if len(reaction_order) <= 1:
        return {}

    return {
        "site_id": site_id,
        "site_hint": _reaction_site_hint(site_candidates),
        "same_site_reaction_count": len(reaction_order),
        "competing_reaction_names": [
            by_reaction[rid] for rid in reaction_order if rid != current_reaction_id
        ],
    }


def explore_reaction(decision_context: Dict[str, Any], reaction_id: str) -> Dict[str, Any]:
    """Return all candidates/templates for one concrete reaction family."""
    units = build_candidate_units(decision_context)
    matching: List[CandidateUnit] = []
    reaction_name = ""
    matched_reaction_id = ""
    for unit in units:
        if unit.source == "terminal":
            continue
        name = _candidate_reaction_name(unit)
        rid = _candidate_reaction_id(unit)
        if rid == reaction_id or name == reaction_id:
            matching.append(unit)
            reaction_name = name
            matched_reaction_id = rid

    if not matching:
        available = [
            {
                "reaction_id": family["reaction_id"],
                "reaction_name": family["reaction_name"],
                "action_count": family["candidate_count"],
            }
            for family in build_reaction_index(decision_context)
        ]
        return {"error": f"unknown reaction_id: {reaction_id}", "available": available}

    sources = _unique_in_order(candidate.source for candidate in matching)
    executable_actions = [candidate for candidate in matching if not _requires_llm_completion(candidate)]
    result = {
        "reaction_id": matched_reaction_id or _reaction_id(reaction_name),
        "reaction_name": reaction_name,
        "sources": sources,
        "action_count": len(matching),
        "actions": [_action_dict(unit, compact=False) for unit in matching],
        "risk_tags": _merge_tags(matching),
        "next_step": "try_action(action_id)" if executable_actions else "propose_action",
    }
    if any(
        "intentional_attachment_placeholder" in candidate.risk_tags
        for candidate in matching
    ):
        result["open_action_requires_completion"] = True
    same_site_context = _same_site_context(
        units,
        matching,
        matched_reaction_id or _reaction_id(reaction_name),
    )
    if same_site_context:
        result["same_site_context"] = same_site_context
    return result


def find_candidate(decision_context: Dict[str, Any], candidate_id: str) -> Optional[Dict[str, Any]]:
    """Find one normalized candidate by stable id."""
    units = build_candidate_units(decision_context)
    for unit in units:
        if unit.candidate_id == candidate_id:
            return unit.to_dict(compact=False)
    for unit in units:
        public_id = _site_aligned_action_id(unit)
        if public_id == candidate_id:
            payload = unit.to_dict(compact=False)
            payload["candidate_id"] = candidate_id
            payload["canonical_candidate_id"] = unit.candidate_id
            return payload
    for unit in units:
        if candidate_id in _legacy_action_ids(unit):
            return unit.to_dict(compact=False)
    return None


def find_action(decision_context: Dict[str, Any], action_id: str) -> Optional[Dict[str, Any]]:
    """Find one normalized action by stable action id."""
    candidate = find_candidate(decision_context, action_id)
    if candidate is None:
        return None
    return _candidate_payload_to_action(candidate)
