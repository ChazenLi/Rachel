"""Prompt-policy and experience-card mounting for Rachel LLM payloads.

This module is intentionally deterministic. It selects short, executable
experience cards from the current molecule/site/sandbox context; it does not
perform chemistry generation and does not expose diagnostic payloads.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional, Set

from Rachel.display_language import (
    DEFAULT_ENGLISH_SELF_PROMPT,
    project_experience_card,
)
from Rachel.chem_tools.cs_score import CS_TRIVIAL
from Rachel.chem_tools.topology_intent import action_declares_topology_change
from Rachel.chem_tools.validation_contract import (
    build_validation_contract,
)
from Rachel.knowledge.conditions import condition_matches


_EXPERIENCE_CARDS_CACHE: Dict[str, List[Dict[str, Any]]] = {}
_PROMPT_POLICY_CACHE: Dict[str, List[Dict[str, Any]]] = {}


def load_experience_cards(knowledge_profile=None) -> List[Dict[str, Any]]:
    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    cached = _EXPERIENCE_CARDS_CACHE.get(knowledge_profile.digest)
    if cached is not None:
        return list(cached)
    data = knowledge_profile.get("prompt.experience_cards")
    cards = data.get("cards", [])
    if not isinstance(cards, list):
        return []
    selected = [card for card in cards if isinstance(card, dict) and card.get("id")]
    _EXPERIENCE_CARDS_CACHE[knowledge_profile.digest] = selected
    return list(selected)


def load_prompt_policy(knowledge_profile=None) -> List[Dict[str, Any]]:
    """Return resolved declarative prompt-policy entries for one profile."""
    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    cached = _PROMPT_POLICY_CACHE.get(knowledge_profile.digest)
    if cached is not None:
        return list(cached)
    data = knowledge_profile.get("prompt.runtime_policy")
    entries = data.get("entries", [])
    if not isinstance(entries, list):
        return []
    selected = [
        entry
        for entry in entries
        if isinstance(entry, dict) and entry.get("id") and entry.get("kind")
    ]
    _PROMPT_POLICY_CACHE[knowledge_profile.digest] = selected
    return list(selected)


def _prompt_policy_entries(kind: str, knowledge_profile=None) -> List[Dict[str, Any]]:
    entries = [
        entry
        for entry in load_prompt_policy(knowledge_profile)
        if entry.get("kind") == kind
    ]
    return sorted(
        entries,
        key=lambda entry: (
            -int(entry.get("priority", 0) or 0),
            str(entry.get("id", "")),
        ),
    )


def _runtime_config(knowledge_profile=None) -> Dict[str, Any]:
    entries = _prompt_policy_entries("runtime_config", knowledge_profile)
    if not entries:
        raise ValueError("prompt runtime policy has no runtime_config entry")
    return entries[0]


def _runtime_config_set(name: str, knowledge_profile=None) -> Set[str]:
    return {
        _normalize_tag(item)
        for item in _runtime_config(knowledge_profile).get(name, []) or []
        if _normalize_tag(item)
    }


def _stage_policy_entry(stage: str, knowledge_profile=None) -> Optional[Dict[str, Any]]:
    matches = [
        entry
        for entry in _prompt_policy_entries("stage_policy", knowledge_profile)
        if str(entry.get("stage", "")) == stage
    ]
    return matches[0] if matches else None


def _event_policy_entries(
    events: List[str], knowledge_profile=None
) -> List[Dict[str, Any]]:
    active = set(events)
    return [
        entry
        for entry in _prompt_policy_entries("event_tags", knowledge_profile)
        if str(entry.get("event", "")) in active
    ]


def _active_policy_entries(
    kind: str,
    *,
    tags: Set[str],
    events: List[str],
    audience: Optional[str] = None,
    knowledge_profile=None,
) -> List[Dict[str, Any]]:
    active_events = set(events)
    return [
        entry
        for entry in _prompt_policy_entries(kind, knowledge_profile)
        if not entry.get("audience") or entry.get("audience") == audience
        if _card_activation_matches(entry, tags, active_events)
    ]


def build_prompt_mount(
    stage: str,
    *,
    decision_context: Optional[Dict[str, Any]] = None,
    payload: Optional[Dict[str, Any]] = None,
    candidate: Optional[Dict[str, Any]] = None,
    attempt: Optional[Dict[str, Any]] = None,
    attempts: Optional[List[Dict[str, Any]]] = None,
    chemist_guidance: Optional[List[Dict[str, Any]]] = None,
    route_plan: Optional[Dict[str, Any]] = None,
    route_strategy: Optional[Dict[str, Any]] = None,
    max_cards: Optional[int] = None,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """Return the short prompt mount for the current LLM interaction stage."""
    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    config = _runtime_config(knowledge_profile)
    stage = _canonical_stage(stage, knowledge_profile)
    prompt_events = derive_prompt_events(
        stage,
        decision_context=decision_context,
        payload=payload,
        candidate=candidate,
        attempt=attempt,
        attempts=attempts,
        chemist_guidance=chemist_guidance,
        route_plan=route_plan,
        route_strategy=route_strategy,
        knowledge_profile=knowledge_profile,
    )
    event_tags = _tags_from_prompt_events(
        prompt_events, knowledge_profile=knowledge_profile
    )
    tags = _extract_tags(
        stage,
        decision_context=decision_context,
        payload=payload,
        candidate=candidate,
        attempt=attempt,
        attempts=attempts,
        chemist_guidance=chemist_guidance,
        route_plan=route_plan,
        route_strategy=route_strategy,
        prompt_events=prompt_events,
        knowledge_profile=knowledge_profile,
    )
    card_limit = _effective_max_cards(
        stage,
        prompt_events=prompt_events,
        tags=tags,
        explicit_max_cards=max_cards,
        knowledge_profile=knowledge_profile,
    )
    active_cards = _select_cards(
        tags,
        event_tags=event_tags,
        prompt_events=prompt_events,
        max_cards=card_limit,
        knowledge_profile=knowledge_profile,
    )
    stage_policy = _stage_policy_entry(stage, knowledge_profile)
    if stage_policy is not None:
        self_prompt = list(stage_policy.get("self_prompts") or [])
        stage_mode = str(stage_policy.get("mode", "audit"))
        command_policy = {
            key: list(value)
            for key, value in (stage_policy.get("command_policy") or {}).items()
            if isinstance(value, list)
        }
    else:
        default_stage = str(config.get("default_stage", "context_compact"))
        default_policy = _stage_policy_entry(default_stage, knowledge_profile) or {}
        self_prompt = list(default_policy.get("self_prompts") or [])
        stage_mode = "audit"
        command_policy = {
            key: list(value)
            for key, value in (config.get("default_command_policy") or {}).items()
            if isinstance(value, list)
        }
    self_prompt_entries = _active_policy_entries(
        "self_prompt",
        tags=tags,
        events=prompt_events,
        audience=stage_mode,
        knowledge_profile=knowledge_profile,
    )
    prepended_prompts = [
        str(entry.get("text", ""))
        for entry in self_prompt_entries
        if entry.get("placement") == "prepend" and str(entry.get("text", "")).strip()
    ]
    appended_prompts = [
        str(entry.get("text", ""))
        for entry in self_prompt_entries
        if entry.get("placement") != "prepend" and str(entry.get("text", "")).strip()
    ]
    self_prompt = list(dict.fromkeys(
        prepended_prompts + self_prompt + appended_prompts
    ))
    audit_field_entries = _active_policy_entries(
        "audit_fields",
        tags=tags,
        events=prompt_events,
        audience=stage_mode,
        knowledge_profile=knowledge_profile,
    )
    allowed_audit_fields = set(config.get("required_audit_fields") or [])
    active_audit_fields = [
        str(field)
        for entry in audit_field_entries
        for field in entry.get("fields", []) or []
        if str(field) in allowed_audit_fields
    ]
    active_audit_fields = list(dict.fromkeys(active_audit_fields))
    tool_hint_entries = _active_policy_entries(
        "tool_hint",
        tags=tags,
        events=prompt_events,
        audience=stage_mode,
        knowledge_profile=knowledge_profile,
    )
    active_tool_hints = [
        str(entry.get("text", ""))
        for entry in tool_hint_entries
        if str(entry.get("text", "")).strip()
    ]
    standing_rule_entries = _prompt_policy_entries(
        "standing_rule", knowledge_profile
    )
    guardrail_entries = _active_policy_entries(
        "quality_guardrail",
        tags=tags,
        events=prompt_events,
        audience=stage_mode,
        knowledge_profile=knowledge_profile,
    )
    event_policy_entries = _event_policy_entries(
        prompt_events, knowledge_profile
    )
    active_quality_guardrails = [
        str(entry.get("text", ""))
        for entry in guardrail_entries
        if str(entry.get("text", "")).strip()
    ]
    active_quality_guardrails = list(dict.fromkeys(active_quality_guardrails))
    guardrail_limit = int(
        (stage_policy or {}).get(
            "guardrail_limit", 5 if stage == "terminal" else 4
        )
        or 0
    )
    active_quality_guardrails = active_quality_guardrails[:guardrail_limit]
    runtime_entries = (
        [config]
        + standing_rule_entries
        + ([stage_policy] if stage_policy is not None else [])
        + event_policy_entries
        + guardrail_entries
        + self_prompt_entries
        + audit_field_entries
        + tool_hint_entries
    )
    policy_refs = [
        knowledge_profile.source("prompt.runtime_policy", str(entry["id"]))
        for entry in runtime_entries
    ]
    card_refs = [
        knowledge_profile.source("prompt.experience_cards", str(card["id"]))
        for card in active_cards
    ]
    brief_policy_refs = []
    for entry in (
        [config]
        + ([stage_policy] if stage_policy is not None else [])
        + event_policy_entries
        + guardrail_entries
        + self_prompt_entries
        + audit_field_entries
        + tool_hint_entries
    ):
        source = knowledge_profile.source(
            "prompt.runtime_policy", str(entry["id"])
        )
        if source.get("pack_kind") != "base":
            brief_policy_refs.append(source)
    return {
        "stage": stage,
        "stage_mode": stage_mode,
        "prompt_state": {
            "stage": stage,
            "events": prompt_events,
        },
        "policy_version": str(config.get("policy_version", "")),
        "card_source": str(config.get("card_source", "")),
        "language_guidance": DEFAULT_ENGLISH_SELF_PROMPT,
        "standing_rules": [
            str(entry.get("text", ""))
            for entry in standing_rule_entries
            if str(entry.get("text", "")).strip()
        ],
        "quality_guardrails": [
            str(entry.get("text", ""))
            for entry in _prompt_policy_entries(
                "quality_guardrail", knowledge_profile
            )
            if str(entry.get("text", "")).strip()
        ],
        "active_quality_guardrails": active_quality_guardrails,
        "command_policy": command_policy,
        "active_experience_card_ids": [card["id"] for card in active_cards],
        "active_experience_cards": [_compact_card(card) for card in active_cards],
        "chemist_guidance": _compact_guidance_list(chemist_guidance),
        "route_plan_brief": _compact_route_plan(route_plan),
        "route_strategy_brief": _compact_route_strategy(route_strategy),
        "self_prompt": self_prompt,
        "required_audit_fields": list(config.get("required_audit_fields") or []),
        "active_audit_fields": active_audit_fields,
        "active_tool_hints": active_tool_hints,
        "matched_tags": sorted(tags),
        "knowledge_profile_hash": knowledge_profile.digest,
        "knowledge_refs": policy_refs + card_refs,
        "brief_knowledge_refs": brief_policy_refs + card_refs,
    }


def build_prompt_brief(mount: Dict[str, Any]) -> Dict[str, Any]:
    """Project an internal prompt mount into the short LLM-facing payload.

    The full mount is useful for deterministic card selection and audit, but it
    repeats policy metadata that should not be read by the model at every step.
    This projection keeps only the stage, event state, next command hints, and
    the compact card text that the model can act on immediately.
    """
    prompt_state = mount.get("prompt_state") or {}
    command_policy = mount.get("command_policy") or {}
    cards = mount.get("active_experience_cards") or []
    stage = str(mount.get("stage", "") or "")
    events = list(prompt_state.get("events") or [])
    self_prompt = list(mount.get("self_prompt") or [])
    if str(mount.get("stage_mode", "")) != "discovery":
        self_prompt = self_prompt[:2]
    brief = {
        "stage": stage,
        "events": events,
        "next_actions": list(command_policy.get("primary_next") or [])[:2],
        "language_guidance": str(mount.get("language_guidance", "")),
        "quality_guardrails": list(mount.get("active_quality_guardrails") or []),
        "active_experience_card_ids": list(mount.get("active_experience_card_ids") or []),
        "experience_prompts": [_compact_card(card) for card in cards],
        "chemist_guidance": list(mount.get("chemist_guidance") or [])[:2],
        "route_plan_brief": _brief_route_plan(mount.get("route_plan_brief") or {}, stage),
        "route_strategy_brief": _brief_route_strategy(
            mount.get("route_strategy_brief") or {},
            stage,
        ),
        "self_prompt": self_prompt,
        "knowledge_profile_hash": str(mount.get("knowledge_profile_hash", "")),
        "knowledge_refs": list(
            mount.get("brief_knowledge_refs")
            or mount.get("knowledge_refs")
            or []
        ),
    }
    audit_fields = list(mount.get("active_audit_fields") or [])
    if audit_fields:
        brief["required_audit_fields"] = audit_fields
    tool_hints = list(mount.get("active_tool_hints") or [])
    if tool_hints:
        brief["tool_hints"] = tool_hints
    return {key: value for key, value in brief.items() if value not in ("", [], {})}


def _brief_route_plan(route_plan: Dict[str, Any], stage: str) -> Dict[str, Any]:
    if not route_plan:
        return {}
    return {
        key: list(value) if isinstance(value, list) else value
        for key, value in route_plan.items()
        if value not in (None, "", [], {})
    }


def _brief_route_strategy(route_strategy: Dict[str, Any], stage: str) -> Dict[str, Any]:
    if not route_strategy or stage not in {
        "context_compact",
        "reaction_sites",
        "route_sketch",
        "propose_action",
        "try_action",
        "sandbox_try",
        "sandbox_list",
        "commit",
        "terminal",
    }:
        return {}
    brief = {
        key: route_strategy.get(key)
        for key in ("id", "macro_strategy", "next_executable_step", "terminal_review")
        if route_strategy.get(key) not in (None, "", [], {}, False)
    }
    if stage in {"route_sketch", "terminal"} and route_strategy.get("problem"):
        brief["problem"] = route_strategy.get("problem")
    rescue_steps = []
    for step in list(route_strategy.get("rescue_steps") or [])[:2]:
        if not isinstance(step, dict):
            continue
        compact_step = {
            key: step.get(key)
            for key in (
                "step_idx",
                "reaction_name",
                "continuation_precursor",
                "status",
            )
            if step.get(key) not in (None, "", [], {})
        }
        if stage in {"route_sketch", "terminal"}:
            for key in ("target_smiles", "target_hint", "expected_precursors"):
                if step.get(key) not in (None, "", [], {}):
                    compact_step[key] = step.get(key)
        if compact_step:
            rescue_steps.append(compact_step)
    if rescue_steps:
        brief["rescue_steps"] = rescue_steps
    return brief


def _canonical_stage(stage: str, knowledge_profile=None) -> str:
    aliases = _runtime_config(knowledge_profile).get("stage_aliases") or {}
    return str(aliases.get(stage, stage))


def _active_card_context_entry(
    tags: Set[str], events: Set[str], knowledge_profile=None
) -> Optional[Dict[str, Any]]:
    contexts = _runtime_config(knowledge_profile).get("selection_contexts") or []
    matches = [
        item
        for item in contexts
        if isinstance(item, dict)
        and item.get("name")
        and condition_matches(
            item.get("activation"), {"tags": tags, "events": events}
        )
    ]
    if not matches:
        return None
    matches.sort(
        key=lambda item: (
            -int(item.get("priority", 0) or 0),
            str(item.get("name", "")),
        )
    )
    return matches[0]


def _active_card_context(
    tags: Set[str], events: Set[str], knowledge_profile=None
) -> str:
    entry = _active_card_context_entry(tags, events, knowledge_profile)
    return str(entry.get("name", "")) if entry else ""


def _effective_max_cards(
    stage: str,
    *,
    prompt_events: List[str],
    tags: Set[str],
    explicit_max_cards: Optional[int],
    knowledge_profile=None,
) -> int:
    if explicit_max_cards is not None:
        return max(0, int(explicit_max_cards))
    active_events = set(prompt_events or [])
    config = _runtime_config(knowledge_profile)
    limits = config.get("card_limits") or {}
    canonical_stage = _canonical_stage(stage, knowledge_profile)
    stage_policy = _stage_policy_entry(canonical_stage, knowledge_profile) or {}
    context_entry = _active_card_context_entry(
        tags, active_events, knowledge_profile
    )
    context = str((context_entry or {}).get("name", ""))
    if (
        stage_policy.get("mode") == "discovery"
        and not (context_entry or {}).get("limit_discovery")
    ):
        return len(load_experience_cards(knowledge_profile))
    if context and context in limits:
        return max(0, int(limits[context]))
    return max(0, int(limits.get("default", 4)))


def _compact_card(card: Dict[str, Any]) -> Dict[str, str]:
    card = project_experience_card(card)
    return {
        "id": str(card.get("id", "")),
        "one_line": str(card.get("one_line", "")),
        "action_prompt": str(card.get("action_prompt", "")),
        "avoid": str(card.get("avoid", "")),
    }


def _compact_guidance_list(guidance: Optional[List[Dict[str, Any]]]) -> List[Dict[str, str]]:
    compact: List[Dict[str, str]] = []
    for item in guidance or []:
        if not isinstance(item, dict):
            continue
        guidance_id = str(item.get("id", "")).strip()
        summary = str(item.get("summary", "")).strip()
        if not guidance_id or not summary:
            continue
        entry = {
            "id": guidance_id,
            "intent": str(item.get("intent", "") or "directive").strip()[:80],
            "summary": summary[:240],
        }
        for key in ("hypothesis_status", "hypothesis_basis"):
            value = str(item.get(key, "")).strip()
            if value:
                entry[key] = value
        compact.append(entry)
        if len(compact) >= 2:
            break
    return compact


def _compact_route_strategy(route_strategy: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not isinstance(route_strategy, dict):
        return {}
    sketch_id = str(route_strategy.get("id", "") or route_strategy.get("route_sketch_id", "")).strip()
    macro = str(route_strategy.get("macro_strategy", "")).strip()
    next_step = str(route_strategy.get("next_executable_step", "")).strip()
    if not sketch_id or not macro:
        return {}
    brief = {
        "id": sketch_id,
        "macro_strategy": macro[:240],
    }
    if next_step:
        brief["next_executable_step"] = next_step[:120]
    problem = str(route_strategy.get("problem", "")).strip()
    if problem:
        brief["problem"] = problem[:160]
    if route_strategy.get("terminal_review"):
        brief["terminal_review"] = True
    rescue_steps = []
    for step in route_strategy.get("rescue_steps", []) or []:
        if not isinstance(step, dict):
            continue
        item = {"step_idx": step.get("step_idx", len(rescue_steps))}
        for key in (
            "target_smiles",
            "target_hint",
            "reaction_name",
            "expected_precursors",
            "continuation_precursor",
            "status",
        ):
            value = step.get(key)
            if value not in ("", [], {}, None):
                item[key] = value
        rescue_steps.append(item)
        if len(rescue_steps) >= 3:
            break
    if rescue_steps:
        brief["rescue_steps"] = rescue_steps
    return brief


def _compact_route_plan(route_plan: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not isinstance(route_plan, dict):
        return {}
    plan_id = str(route_plan.get("id", "")).strip()
    thesis = str(route_plan.get("route_thesis", "")).strip()
    if not plan_id or not thesis:
        return {}
    brief: Dict[str, Any] = {
        "id": plan_id,
        "revision": int(route_plan.get("revision", 0) or 0),
        "route_thesis": thesis,
    }
    route_mode = str(route_plan.get("route_mode", "")).strip()
    if route_mode:
        brief["route_mode"] = route_mode
    for key in (
        "key_disconnections",
        "preferred_precursor_logic",
        "protect_or_preserve",
        "mode_evidence",
        "strategic_risks",
        "revision_triggers",
    ):
        raw_items = route_plan.get(key, []) or []
        if not isinstance(raw_items, list):
            continue
        items = [
            str(item).strip()
            for item in raw_items
            if str(item).strip()
        ]
        if items:
            brief[key] = items
    terminal_policy = str(route_plan.get("terminal_rescue_policy", "")).strip()
    if terminal_policy:
        brief["terminal_rescue_policy"] = terminal_policy
    reason = str(route_plan.get("last_revision_reason", "")).strip()
    if reason:
        brief["last_revision_reason"] = reason
    return brief


def derive_prompt_events(
    stage: str,
    *,
    decision_context: Optional[Dict[str, Any]] = None,
    payload: Optional[Dict[str, Any]] = None,
    candidate: Optional[Dict[str, Any]] = None,
    attempt: Optional[Dict[str, Any]] = None,
    attempts: Optional[List[Dict[str, Any]]] = None,
    chemist_guidance: Optional[List[Dict[str, Any]]] = None,
    route_plan: Optional[Dict[str, Any]] = None,
    route_strategy: Optional[Dict[str, Any]] = None,
    knowledge_profile=None,
) -> List[str]:
    """Derive compact prompt-state events from structured runtime payloads."""
    stage = _canonical_stage(stage, knowledge_profile)
    events: List[str] = []
    seen: Set[str] = set()

    _add_event(events, seen, f"stage.{_event_suffix(stage)}")
    if stage == "commit":
        _add_event(events, seen, "decision.commit_requested")

    _add_events_from_payload(events, seen, payload)
    _add_events_from_candidate(events, seen, candidate)
    _add_events_from_attempt(events, seen, attempt)
    if stage != "commit":
        for item in attempts or []:
            _add_events_from_attempt(events, seen, item)
    _add_events_from_chemist_guidance(events, seen, chemist_guidance)
    _add_events_from_route_plan(events, seen, route_plan, stage=stage)
    _add_events_from_route_strategy(events, seen, route_strategy)
    if stage != "commit":
        _add_events_from_attempt_set(events, seen, attempts)
    _add_events_from_decision_context(
        events,
        seen,
        decision_context,
        stage=stage,
        route_plan_active=bool(_compact_route_plan(route_plan)),
    )

    return events


def _select_cards(
    tags: Set[str],
    *,
    event_tags: Optional[Set[str]] = None,
    prompt_events: Optional[List[str]] = None,
    max_cards: int,
    knowledge_profile=None,
) -> List[Dict[str, Any]]:
    scored: List[tuple[int, int, Dict[str, Any]]] = []
    active_events = set(prompt_events or [])
    event_tag_set = event_tags or set()
    facts = {"tags": tags, "events": active_events}
    for order, card in enumerate(load_experience_cards(knowledge_profile)):
        if not _card_activation_matches(card, tags, active_events):
            continue
        selection = card.get("selection") or {}
        if not isinstance(selection, dict):
            selection = {}
        if not condition_matches(selection.get("eligibility"), facts):
            continue
        card_tags = {_normalize_tag(tag) for tag in card.get("tags", [])}
        overlap = card_tags & tags
        event_overlap = card_tags & event_tag_set
        trigger_overlap = {
            str(trigger)
            for trigger in card.get("triggers", []) or []
            if str(trigger) in active_events
        }
        if not overlap and not trigger_overlap:
            continue
        score = (
            len(overlap) * 10
            + len(event_overlap) * 5
            + len(trigger_overlap) * 15
            + int(card.get("priority", 0) or 0)
        )
        for boost in selection.get("boosts", []) or []:
            if not isinstance(boost, dict):
                continue
            if condition_matches(boost.get("when"), facts):
                score += int(boost.get("score", 0) or 0)
        scored.append((score, -order, card))
    scored.sort(reverse=True)
    ranked = [card for _, _, card in scored]
    context = _active_card_context(tags, active_events, knowledge_profile)
    if context:
        return _select_context_cards(
            ranked,
            context=context,
            facts=facts,
            max_cards=max_cards,
        )
    return ranked[:max_cards]


def _card_activation_matches(
    card: Dict[str, Any],
    tags: Set[str],
    active_events: Set[str],
) -> bool:
    activation = card.get("activation") or {}
    return condition_matches(
        activation,
        {
            "tags": tags,
            "events": active_events,
        },
    )


def _has_explicit_topology_signal(tags: Set[str]) -> bool:
    return "topology_signal" in tags


def _select_context_cards(
    ranked: List[Dict[str, Any]],
    *,
    context: str,
    facts: Dict[str, Any],
    max_cards: int,
) -> List[Dict[str, Any]]:
    if max_cards <= 0:
        return []
    selected: List[Dict[str, Any]] = []
    selected_ids: Set[str] = set()
    slotted: List[tuple[int, int, Dict[str, Any]]] = []
    for rank, card in enumerate(ranked):
        matching_orders = [
            int(slot.get("order", 0) or 0)
            for slot in (card.get("selection") or {}).get("slots", []) or []
            if isinstance(slot, dict)
            and str(slot.get("context", "")) == context
            and condition_matches(slot.get("when"), facts)
        ]
        if matching_orders:
            slotted.append((min(matching_orders), rank, card))

    for _, _, card in sorted(slotted, key=lambda item: (item[0], item[1])):
        card_id = str(card.get("id", ""))
        if card_id and card_id not in selected_ids:
            selected.append(card)
            selected_ids.add(card_id)
        if len(selected) >= max_cards:
            return selected[:max_cards]

    for card in ranked:
        card_id = str(card.get("id", ""))
        if card_id in selected_ids:
            continue
        selected.append(card)
        selected_ids.add(card_id)
        if len(selected) >= max_cards:
            break
    return selected


def _extract_tags(
    stage: str,
    *,
    decision_context: Optional[Dict[str, Any]],
    payload: Optional[Dict[str, Any]],
    candidate: Optional[Dict[str, Any]],
    attempt: Optional[Dict[str, Any]],
    attempts: Optional[List[Dict[str, Any]]],
    chemist_guidance: Optional[List[Dict[str, Any]]] = None,
    route_plan: Optional[Dict[str, Any]] = None,
    route_strategy: Optional[Dict[str, Any]] = None,
    prompt_events: Optional[List[str]] = None,
    knowledge_profile=None,
) -> Set[str]:
    tags = {_normalize_tag(stage)}
    tags.update(
        _tags_from_prompt_events(
            prompt_events or [], knowledge_profile=knowledge_profile
        )
    )

    text_parts: List[str] = []
    text_parts.extend(_collect_molecule_context_text(decision_context))
    text_parts.extend(_collect_local_action_text(payload))
    text_parts.extend(_collect_local_action_text(candidate))
    text_parts.extend(_collect_local_action_text(attempt))
    text_parts.extend(_collect_text(_compact_guidance_list(chemist_guidance)))

    text = " ".join(text_parts).lower()
    tags.update(_derive_tags_from_text(text))

    if payload:
        tags.update(_derive_tags_from_payload(payload, knowledge_profile))
    if candidate:
        tags.update(_derive_tags_from_candidate(candidate, knowledge_profile))
    if attempt:
        tags.update(_derive_tags_from_attempt(attempt, knowledge_profile))
    if stage != "commit":
        for item in attempts or []:
            tags.update(_derive_tags_from_attempt(item, knowledge_profile))
    if _compact_route_plan(route_plan):
        tags.add("route_plan")
        route_mode = _normalize_tag(route_plan.get("route_mode", "")) if route_plan else ""
        if route_mode:
            tags.update({"route_mode", route_mode})
    if _compact_route_strategy(route_strategy):
        tags.update({"route_sketch", "strategic_rescue"})

    return {tag for tag in tags if tag}


def _tags_from_prompt_events(
    prompt_events: List[str], *, knowledge_profile=None
) -> Set[str]:
    tags: Set[str] = set()
    for event in prompt_events:
        tags.add(_normalize_tag(event))
        tags.add(_event_suffix(event))
    for entry in _event_policy_entries(prompt_events, knowledge_profile):
        tags.update(
            _normalize_tag(tag)
            for tag in entry.get("tags", []) or []
            if _normalize_tag(tag)
        )
    return {tag for tag in tags if tag}


def _add_events_from_chemist_guidance(
    events: List[str],
    seen: Set[str],
    guidance: Optional[List[Dict[str, Any]]],
) -> None:
    active = [item for item in (guidance or []) if isinstance(item, dict)]
    if not active:
        return
    _add_event(events, seen, "chemist.directive")
    for item in active:
        intent = str(item.get("intent", "")).lower()
        if item.get("site_hint"):
            _add_event(events, seen, "chemist.site_hint")
        if item.get("reaction_hint"):
            _add_event(events, seen, "chemist.reaction_hint")
        if item.get("precursors"):
            _add_event(events, seen, "chemist.precursor_hint")
        if item.get("terminal_hint") or "terminal" in intent:
            _add_event(events, seen, "chemist.terminal_hint")
        if "approval" in intent or "approve" in intent:
            _add_event(events, seen, "chemist.approval")


def _add_events_from_route_strategy(
    events: List[str],
    seen: Set[str],
    route_strategy: Optional[Dict[str, Any]],
) -> None:
    if not isinstance(route_strategy, dict) or not _compact_route_strategy(route_strategy):
        return
    _add_event(events, seen, "strategy.route_sketch_active")
    if route_strategy.get("terminal_review"):
        _add_event(events, seen, "strategy.advanced_terminal_rescue_requested")
    next_step = str(route_strategy.get("next_executable_step", "")).lower()
    if "propose" in next_step or "custom" in next_step:
        _add_event(events, seen, "strategy.route_sketch_requested")


def _add_events_from_route_plan(
    events: List[str],
    seen: Set[str],
    route_plan: Optional[Dict[str, Any]],
    *,
    stage: str,
) -> None:
    if not isinstance(route_plan, dict) or not _compact_route_plan(route_plan):
        return
    _add_event(events, seen, "strategy.route_plan_active")
    if stage == "route_plan" and str(route_plan.get("route_mode", "") or "").strip():
        _add_event(events, seen, "strategy.route_mode_triage")
    reason = str(route_plan.get("last_revision_reason", "")).lower()
    try:
        revision = int(route_plan.get("revision", 0) or 0)
    except (TypeError, ValueError):
        revision = 0
    if stage == "route_plan" and (revision > 0 or reason not in {"", "initial"}):
        _add_event(events, seen, "strategy.route_plan_revised")


def _add_events_from_attempt_set(
    events: List[str],
    seen: Set[str],
    attempts: Optional[List[Dict[str, Any]]],
) -> None:
    if not attempts:
        return
    weak_count = 0
    for attempt in attempts:
        validation = build_validation_contract(
            attempt.get("forward_validation") or {},
            validation_micro=attempt.get("validation_micro") or {},
            site_audit=attempt.get("site_audit") or {},
            eas_site_audit=attempt.get("eas_site_audit") or {},
            execution_success=attempt.get("success"),
        )
        state = str((validation.get("decision_gate") or {}).get("state", ""))
        if attempt.get("success") is False or state in {"blocked", "proof_required"}:
            weak_count += 1
    if weak_count >= 2:
        _add_event(events, seen, "strategy.route_sketch_requested")


_LOCAL_ACTION_TEXT_KEYS = {
    "source",
    "source_label",
    "reaction_type",
    "reaction_name",
    "action_label",
    "candidate_label",
    "site_type",
    "site_hint",
    "risk_hint",
    "risk_tags",
    "role_pair",
    "bond_type",
    "bond_fg_context",
    "rationale_summary",
    "why_existing_actions_rejected",
    "why_existing_candidates_rejected",
    "intended_deltas",
    "expected_ring_change",
    "mechanistic_evidence",
    "reaction",
}


def _collect_local_action_text(value: Any) -> Iterable[str]:
    if not isinstance(value, dict):
        return []
    summary = value.get("action_summary") or value.get("candidate_summary") or {}
    source = summary if isinstance(summary, dict) and summary else value
    texts: List[str] = []
    for key in _LOCAL_ACTION_TEXT_KEYS:
        if key in source:
            texts.extend(_collect_text(source.get(key)))
    return texts


def _collect_molecule_context_text(value: Any) -> Iterable[str]:
    if not isinstance(value, dict):
        return []
    texts: List[str] = []
    for key in (
        "smiles",
        "molecule",
        "scaffold",
        "scaffold_tags",
        "functional_groups",
        "functional_group_brief",
        "fg_summary",
    ):
        if key in value:
            texts.extend(_collect_text(value.get(key)))
    return texts


def _collect_text(value: Any) -> Iterable[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        texts: List[str] = []
        for item in value.values():
            texts.extend(_collect_text(item))
        return texts
    if isinstance(value, list):
        texts = []
        for item in value:
            texts.extend(_collect_text(item))
        return texts
    return []


def _precursor_normalization_tags(raw: Any) -> Set[str]:
    if not isinstance(raw, dict):
        return set()
    if not raw.get("changed") and not raw.get("normalizations"):
        return set()
    return {
        "precursor_normalization",
        "organometallic_source_obligation",
        "organometallic",
        "metal_source",
    }


def _derive_tags_from_text(text: str) -> Set[str]:
    tags: Set[str] = set()
    if "electron-poor" in text or "electron_poor" in text or "electron poor" in text:
        tags.add("electron_poor")
    if "pyridine" in text or "pyridyl" in text or "nicotinate" in text or "pyridone" in text:
        tags.update({"pyridine", "heteroaryl"})
    if "pyridone" in text or "lactam" in text or "electronic-state" in text or "electronic state" in text:
        tags.update({"electronic_state_strategy", "scaffold_assembly", "route_mode"})
    if "scaffold assembly" in text or "ring construction" in text:
        tags.update({"scaffold_assembly", "route_mode"})
    if any(
        term in text
        for term in (
            "ring_closure",
            "ring closure",
            "ring opening",
            "ring_opening",
            "scaffold_edit",
            "scaffold edit",
            "new_fused_ring_system",
            "new_spiro_center",
            "new_bridged_system",
            "new_medium_ring",
            "new_macrocycle",
            "high_risk_topology_requires_independent_evidence",
        )
    ):
        tags.update({"topology", "ring_construction", "scaffold_edit"})
    if "fused" in text or "fused_ring" in text:
        tags.add("fused_scaffold")
        if any(term in text for term in ("heteroaryl", "pyridine", "aza", "indole")):
            tags.add("fused_heteroaryl")
    if "spiro" in text:
        tags.add("spiro_scaffold")
    if "bridged" in text:
        tags.add("bridged_scaffold")
    if "macrocycle" in text:
        tags.add("macrocycle_scaffold")
    if "late functionalization" in text or "late fgi" in text or "functional group interconversion" in text:
        tags.update({"late_functionalization", "route_mode"})
    if "terminal hides core" in text or "hidden core" in text or "core construction" in text:
        tags.update({"terminal_hides_core_problem", "advanced_terminal", "route_plan"})
    if "snar" in text or "n-nucleophile" in text:
        tags.update({"snar", "n_nucleophile", "heteroaryl"})
    if "paal" in text and "knorr" in text:
        tags.update({"paal-knorr", "ring_construction", "template_risk", "advanced_terminal", "deep_disconnection"})
    if "suzuki" in text:
        tags.update({"suzuki", "sp2_sp2", "convergent"})
    if "miyaura" in text or "borylation" in text:
        tags.update({"handle_timing", "suzuki"})
    if "halogenation" in text or "aryl_halide" in text or "fluoride" in text or "fluoro" in text:
        tags.update({"handle_timing", "halogenation"})
    if "heteroaryl_fluoride" in text or "heteroaryl fluoride" in text:
        tags.add("heteroaryl_fluoride")
    if (
        "pyridine" in tags
        and any(
            term in text
            for term in (
                "electron poor",
                "electron_poor",
                "chloro",
                "chloride",
                "nitro",
                "cyano",
                "ester",
                "difluoromethyl",
                "trifluoromethyl",
                "halogenation",
            )
        )
    ):
        tags.update({"electron_poor", "electron_poor_pyridine", "heteroaryl"})
    if "heteroaryl" in text or "azaindole" in text or "indole" in text or "pyrrole" in text:
        tags.update({"heteroaryl", "fused_heteroaryl"})
    if "amide" in text:
        tags.update({"amide", "atom_accounting", "convergent"})
    if "protection" in text or "boc" in text or "free amine" in text:
        tags.update({"protection", "compatibility", "free_amine"})
    if "skeleton_imbalance" in text:
        tags.update({"skeleton_imbalance", "multicomponent"})
    if (
        "custom_precursors" in text
        or "llm_proposed" in text
        or "propose_action" in text
        or "propose_candidate" in text
    ):
        tags.update({"custom_precursor", "llm_proposed", "action_rejection"})
    if "proof_required" in text or "proof obligation" in text or "proof_obligation" in text:
        tags.update({"forward_validation", "proof_required", "proof_obligation"})
    if (
        "precursor_normalization" in text
        or "organometallic_source_obligation" in text
        or "organometallic_precursor_normalized" in text
    ):
        tags.update({
            "precursor_normalization",
            "organometallic_source_obligation",
            "organometallic",
        })
    if "atom_mapping" in text:
        tags.update({"atom_mapping", "atom_source"})
    if "formed_bonds" in text:
        tags.update({"atom_mapping", "formed_bonds", "atom_source"})
    if "cleaved_bonds" in text:
        tags.update({"atom_mapping", "cleaved_bonds", "atom_source"})
    if "unmapped_product_atoms" in text:
        tags.update({"atom_mapping", "unmapped_product_atoms", "atom_source"})
    return tags


def _derive_tags_from_payload(
    payload: Dict[str, Any], knowledge_profile=None
) -> Set[str]:
    tags: Set[str] = set()
    for action in payload.get("actions", []) or []:
        if not isinstance(action, dict):
            continue
        local_text = " ".join(_collect_local_action_text(action)).lower()
        tags.update(_derive_tags_from_text(local_text))
        tags.update(_derive_tags_from_candidate(action, knowledge_profile))
    risk_tags = {
        _normalize_tag(tag)
        for tag in payload.get("risk_tags", []) or []
        if _normalize_tag(tag)
    }
    tags.update(risk_tags)
    if _runtime_config_set("topology_signal_risk_tags", knowledge_profile) & risk_tags:
        tags.add("topology_signal")
    if payload.get("in_ring"):
        tags.update({"topology_signal", "topology", "ring_bond", "ring_construction"})
    if payload.get("terminal_review") is True:
        tags.update({"advanced_terminal", "terminal_rescue", "route_sketch"})
    return tags


def _derive_tags_from_candidate(
    candidate: Dict[str, Any], knowledge_profile=None
) -> Set[str]:
    tags: Set[str] = set()
    source = str(candidate.get("source", ""))
    if source == "custom_precursors":
        tags.update({"custom_precursor", "llm_proposed", "action_rejection", "audit"})
    if candidate.get("route_sketch_id") or candidate.get("route_strategy_brief"):
        tags.update({"route_sketch", "strategic_rescue"})
    if action_declares_topology_change(candidate):
        tags.update({"topology_signal", "topology", "ring_construction", "scaffold_edit"})
        expected = str(candidate.get("expected_ring_change", "") or "").lower()
        if "spiro" in expected:
            tags.add("spiro")
        if "fused" in expected:
            tags.add("fused_ring")
        if "bridged" in expected:
            tags.add("bridged")
        if "macro" in expected:
            tags.add("macrocycle")
        for delta in candidate.get("intended_deltas", []) or []:
            normalized = _normalize_tag(delta)
            if normalized:
                tags.add(normalized)
    tags.update(_precursor_normalization_tags(candidate.get("precursor_normalization")))
    risk_tags = {
        _normalize_tag(tag)
        for tag in candidate.get("risk_tags", []) or []
        if _normalize_tag(tag)
    }
    tags.update(risk_tags)
    if _runtime_config_set("topology_signal_risk_tags", knowledge_profile) & risk_tags:
        tags.add("topology_signal")
    return tags


def _derive_tags_from_attempt(
    attempt: Dict[str, Any], knowledge_profile=None
) -> Set[str]:
    tags = {"sandbox"}
    validation = build_validation_contract(
        attempt.get("forward_validation") or {},
        validation_micro=attempt.get("validation_micro") or {},
        site_audit=attempt.get("site_audit") or {},
        eas_site_audit=attempt.get("eas_site_audit") or {},
        execution_success=attempt.get("success"),
    )
    decision_gate = validation.get("decision_gate") or {}
    state = _normalize_tag(decision_gate.get("state", ""))
    if state and state != "unknown":
        tags.update({"forward_validation", state})
    if state == "blocked":
        tags.add("hard_fail")
    if state == "proof_required":
        tags.add("proof_obligation")
    for bucket_name, bucket_tag in (
        ("contradictions", "hard_fail"),
        ("proof_obligations", "proof_obligation"),
        ("evidence_gaps", "evidence_gap"),
        ("tool_limits", "tool_limit"),
        ("warnings", "warning"),
        ("system_errors", "system_error"),
    ):
        for item in validation.get(bucket_name, []) or []:
            code = _normalize_tag(item.get("code", ""))
            if code:
                tags.update({code, bucket_tag})
    summary = attempt.get("action_summary") or attempt.get("candidate_summary") or {}
    if summary.get("source") == "custom_precursors":
        tags.update({"custom_precursor", "llm_proposed", "action_rejection"})
    if isinstance(summary, dict):
        tags.update(_derive_tags_from_candidate(summary, knowledge_profile))
    tags.update(_precursor_normalization_tags(attempt.get("precursor_normalization")))
    tags.update(_precursor_normalization_tags(summary.get("precursor_normalization")))
    observations = validation.get("observations") or {}
    ring_deltas = list(observations.get("ring_deltas", []) or [])
    if ring_deltas:
        tags.update({"topology_signal", "topology", "ring_construction", "scaffold_edit"})
        tags.update(_normalize_tag(code) for code in ring_deltas if _normalize_tag(code))
    if observations.get("atom_mapping"):
        tags.update({"atom_mapping", "atom_source"})
    mechanism = validation.get("mechanism_interpretation") or {}
    if mechanism.get("label") and mechanism.get("state") != "unregistered_family":
        tags.update({"mechanism_interpretation", "reaction_topology"})
    return tags


def _add_event(events: List[str], seen: Set[str], event: str) -> None:
    event = event.strip().lower()
    if event and event not in seen:
        events.append(event)
        seen.add(event)


def _add_events_from_decision_context(
    events: List[str],
    seen: Set[str],
    decision_context: Optional[Dict[str, Any]],
    *,
    stage: str,
    route_plan_active: bool,
) -> None:
    if not decision_context:
        return
    text = " ".join(_collect_molecule_context_text(decision_context)).lower()
    if (
        stage in {"context_compact", "reaction_sites"}
        and not route_plan_active
        and
        any(term in text for term in ("pyridine", "pyridyl", "nicotinate"))
        and any(term in text for term in ("electron poor", "electron_poor", "chloro", "chloride", "ester", "nitro", "cyano", "difluoromethyl"))
    ):
        _add_event(events, seen, "strategy.route_mode_triage")
    if stage in {"reaction_sites", "explore_site"} and (
        "heteroaryl" in text or "azaindole" in text or "indole" in text
    ):
        _add_event(events, seen, "site.fused_heteroaryl_site_sensitive")
    action_count = 0
    for bond in decision_context.get("disconnectable_bonds", []) or []:
        action_count += len(bond.get("alternatives", []) or [])
        action_count += len(bond.get("smart_capping", []) or [])
    action_count += len(decision_context.get("fgi_options", []) or [])
    action_count += len(decision_context.get("custom_candidates", []) or [])
    complexity = decision_context.get("complexity") or {}
    try:
        cs_score = float(complexity.get("cs_score", 0) or 0)
    except (TypeError, ValueError):
        cs_score = 0.0
    if stage in {"context_compact", "reaction_sites"} and action_count <= 1 and cs_score >= CS_TRIVIAL:
        _add_event(events, seen, "strategy.action_space_weak")
        _add_event(events, seen, "strategy.route_sketch_requested")


def _add_events_from_payload(
    events: List[str],
    seen: Set[str],
    payload: Optional[Dict[str, Any]],
) -> None:
    if not payload:
        return
    if payload.get("validation_override"):
        _add_event(events, seen, "decision.commit_with_override")
    if payload.get("terminal_review"):
        _add_event(events, seen, "strategy.advanced_terminal_rescue_requested")
    if payload.get("site_audit"):
        _add_events_from_site_audit(events, seen, payload.get("site_audit") or {})
    site_map = payload.get("site_reaction_map") or []
    for site in site_map if isinstance(site_map, list) else []:
        if (site.get("reaction_count") or 0) > 1 or site.get("competition_hint"):
            _add_event(events, seen, "site.same_site_competition")
    if (payload.get("reaction_count") or 0) > 1 or payload.get("competition_hint"):
        _add_event(events, seen, "site.same_site_competition")
    if payload.get("in_ring"):
        _add_event(events, seen, "site.ring_bond_review")
    hint_text = " ".join(_collect_local_action_text(payload)).lower()
    if payload.get("same_core") or payload.get("site_retentive"):
        _add_event(events, seen, "site.same_core_custom_precursor")
    if payload.get("site_anchor_drift"):
        _add_event(events, seen, "site.site_anchor_drift")
    if payload.get("action_space_weak") or payload.get("route_incoherent"):
        _add_event(events, seen, "strategy.route_sketch_requested")
    if payload.get("open_action_requires_completion"):
        _add_event(events, seen, "action.intentional_attachment_placeholder")
        _add_event(events, seen, "strategy.route_sketch_requested")


def _add_events_from_candidate(
    events: List[str],
    seen: Set[str],
    candidate: Optional[Dict[str, Any]],
) -> None:
    if not candidate:
        return
    summary = candidate.get("action_summary") or candidate.get("candidate_summary") or {}
    source = str(candidate.get("source") or summary.get("source") or "")
    action_id = str(candidate.get("action_id") or candidate.get("candidate_id") or "")
    if source in {"bond", "fgi"}:
        _add_event(events, seen, "action.system_template")
    if source == "smart_capping":
        _add_event(events, seen, "action.smart_capping")
    if source == "terminal":
        _add_event(events, seen, "action.terminal_acceptance")
    if source == "custom_precursors" or action_id.startswith("custom:"):
        _add_event(events, seen, "action.custom_precursors")
    if _precursor_normalization_tags(candidate.get("precursor_normalization") or summary.get("precursor_normalization")):
        _add_event(events, seen, "action.precursor_normalization")
    if candidate.get("route_sketch_id") or summary.get("route_sketch_id"):
        _add_event(events, seen, "strategy.route_sketch_used_for_custom_action")
    risk_tags = {
        _normalize_tag(tag)
        for tag in candidate.get("risk_tags", []) or summary.get("risk_tags", []) or []
    }
    if "intentional_attachment_placeholder" in risk_tags:
        _add_event(events, seen, "action.intentional_attachment_placeholder")
    if (
        candidate.get("duplicate_of")
        or summary.get("duplicate_of")
        or "duplicates_existing_candidate" in risk_tags
        or "duplicates_existing_action" in risk_tags
    ):
        _add_event(events, seen, "action.duplicates_existing_action")


def _add_events_from_attempt(
    events: List[str],
    seen: Set[str],
    attempt: Optional[Dict[str, Any]],
) -> None:
    if not attempt:
        return
    if "success" in attempt:
        _add_event(
            events,
            seen,
            "sandbox.success" if attempt.get("success") else "sandbox.failure",
        )
    _add_events_from_candidate(
        events,
        seen,
        attempt.get("action_summary") or attempt.get("candidate_summary") or attempt,
    )
    if _precursor_normalization_tags(attempt.get("precursor_normalization")):
        _add_event(events, seen, "action.precursor_normalization")
    if attempt.get("site_audit"):
        _add_events_from_site_audit(events, seen, attempt.get("site_audit") or {})

    validation = build_validation_contract(
        attempt.get("forward_validation") or {},
        validation_micro=attempt.get("validation_micro") or {},
        site_audit=attempt.get("site_audit") or {},
        eas_site_audit=attempt.get("eas_site_audit") or {},
        execution_success=attempt.get("success"),
    )
    state = _event_suffix((validation.get("decision_gate") or {}).get("state", ""))
    if state and state != "unknown":
        _add_event(events, seen, f"validation.{state}")
    for bucket_name in (
        "contradictions",
        "proof_obligations",
        "evidence_gaps",
        "tool_limits",
        "warnings",
        "system_errors",
    ):
        for item in validation.get(bucket_name, []) or []:
            code = _event_suffix(item.get("code", ""))
            if code:
                _add_event(events, seen, f"validation.{code}")


def _add_events_from_site_audit(
    events: List[str],
    seen: Set[str],
    site_audit: Dict[str, Any],
) -> None:
    summary = str(site_audit.get("summary", "") or "").lower()
    if "drift" in summary or site_audit.get("site_anchor_drift"):
        _add_event(events, seen, "site.site_anchor_drift")
    if "missing" in summary or site_audit.get("missing_evidence"):
        _add_event(events, seen, "site.site_retention_missing_evidence")


def _event_suffix(value: Any) -> str:
    text = _normalize_tag(value)
    return text.replace("_", "_")


def _normalize_tag(value: Any) -> str:
    text = str(value or "").strip().lower()
    text = text.replace("-", "_")
    text = re.sub(r"[^a-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    if text == "paal_knorr":
        return "paal-knorr"
    return text
