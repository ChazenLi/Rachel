"""Deterministic packet and audit helpers for current-node scouting."""

from __future__ import annotations

import copy
import hashlib
import json
import re
from typing import Any, Dict, List


SCHEMA_VERSION = "rachel.node_scouting.v1"
MIN_TASKS = 5
MAX_TASKS = 10
METHODOLOGY_LENSES = {
    "photo_electrocatalysis": (
        "Test photochemical, electrochemical, or photoelectrochemical activation, "
        "including SET logic and a closed redox or catalytic cycle."
    ),
    "biocatalysis": (
        "Test enzyme, whole-cell, or cofactor-mediated chemistry for useful chemo-, "
        "regio-, or stereocontrol."
    ),
    "metal_dual_catalysis": (
        "Test transition-metal or dual-catalytic bond construction, including "
        "cross-coupling, C-H activation, and cooperative activation."
    ),
    "asymmetric_catalysis": (
        "Test catalytic stereocontrol for enantioselective or diastereoselective "
        "bond construction and identify the source of selectivity."
    ),
    "radical_polarity_special_activation": (
        "Test radical, polarity-reversal, radical-polar crossover, HAT, or other "
        "special activation modes with explicit selectivity and termination logic."
    ),
    "skeletal_editing_rearrangement_cascade": (
        "Test skeletal editing, atom insertion/deletion/exchange, rearrangement, "
        "ring editing, or a controlled cascade at the current structure."
    ),
}
HYPOTHESIS_STATUSES = {
    "precedent_supported",
    "analogy_based",
    "unprecedented_hypothesis",
}
FRONTIER_HYPOTHESIS_STATUS = "unprecedented_hypothesis"
FRONTIER_ANALYSIS_CONTRACT = {
    "indexed_transformation": (
        "Name the proposed bond disconnection or FGI on packet atom/bond indices "
        "and the reasonable direct-precursor relationship."
    ),
    "mechanism_and_balance": (
        "Give the stepwise mechanistic hypothesis, atom/electron/proton balance, "
        "and catalyst or mediator turnover/regeneration."
    ),
    "selectivity": (
        "Explain chemo-, regio-, and stereocontrol, including the proposed source "
        "of asymmetric induction when relevant."
    ),
    "proposed_operation": (
        "Propose testable reagents/catalysts, order of addition, medium, temperature, "
        "atmosphere, light or potential as applicable; label unsupported details."
    ),
    "workup_and_purification": (
        "Propose a compatible quench, workup, and purification strategy, including "
        "likely instability or separation risks."
    ),
    "failure_and_falsification": (
        "Name likely failure modes, diagnostic observations, decisive controls, and "
        "what result would falsify or narrow the hypothesis."
    ),
}
SCOUTING_PROMPT_GUIDANCE = (
    "Use the latest scouting review as advisory current-node hypotheses. "
    "Compare distinct reaction ideas and their stated target sites, collapse "
    "duplicates to one direction, and independently adopt, revise, or ignore "
    "each through the normal Rachel action and sandbox flow. Treat "
    "any unprecedented_hypothesis as explicitly speculative and independently "
    "review its mechanism and practical test before realization. Treat "
    "scouting_source as provenance only; validation comes from the realized action."
)


def _normalized_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _digest(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def compact_action_space(site_result: Dict[str, Any]) -> Dict[str, Any]:
    """Project the first-layer menu without prompt or mutation affordances."""
    rows: List[Dict[str, Any]] = []
    for raw_site in site_result.get("site_reaction_map", []) or []:
        if not isinstance(raw_site, dict):
            continue
        site = {
            key: copy.deepcopy(raw_site.get(key))
            for key in (
                "site_id",
                "site_type",
                "site_hint",
                "action_count",
                "reaction_count",
                "competition_hint",
                "risk_hint",
            )
            if raw_site.get(key) not in (None, "", [], {})
        }
        reactions = []
        for raw_reaction in raw_site.get("reactions", []) or []:
            if not isinstance(raw_reaction, dict):
                continue
            reactions.append({
                key: copy.deepcopy(raw_reaction.get(key))
                for key in (
                    "reaction_id",
                    "reaction_name",
                    "action_count",
                    "source_summary",
                    "risk_hint",
                )
                if raw_reaction.get(key) not in (None, "", [], {})
            })
        site["reactions"] = reactions
        rows.append(site)
    return {"site_reaction_map": rows}


def _normalize_tasks(
    tasks: Any,
    *,
    expansion_reason: str,
    frontier_reason: str,
) -> List[Dict[str, Any]]:
    if not isinstance(tasks, list):
        raise ValueError("tasks must be a list")
    if not MIN_TASKS <= len(tasks) <= MAX_TASKS:
        raise ValueError("tasks must contain 5 to 10 entries")

    normalized: List[Dict[str, Any]] = []
    seen_focuses = set()
    blind_count = 0
    informed_count = 0
    frontier_count = 0
    for index, task in enumerate(tasks):
        if not isinstance(task, dict):
            raise ValueError(f"tasks[{index}] must be an object")
        visibility = _normalized_text(task.get("visibility")).lower()
        if visibility not in {"blind", "informed"}:
            raise ValueError(f"tasks[{index}].visibility must be blind or informed")
        focus = _normalized_text(task.get("focus"))
        if not focus:
            raise ValueError(f"tasks[{index}].focus required")
        focus_key = focus.casefold()
        if focus_key in seen_focuses:
            raise ValueError("task focus values must be unique after normalization")
        seen_focuses.add(focus_key)
        methodology_lens = ""
        if "methodology_lens" in task:
            methodology_lens = _normalized_text(task.get("methodology_lens")).lower()
            if methodology_lens not in METHODOLOGY_LENSES:
                raise ValueError(
                    f"tasks[{index}].methodology_lens must be one of the six supported lenses"
                )
        frontier_hypothesis = task.get("frontier_hypothesis", False)
        if "frontier_hypothesis" in task and not isinstance(frontier_hypothesis, bool):
            raise ValueError(f"tasks[{index}].frontier_hypothesis must be boolean")
        if frontier_hypothesis and not methodology_lens:
            raise ValueError(
                f"tasks[{index}].methodology_lens required for a frontier hypothesis"
            )

        blind_count += visibility == "blind"
        informed_count += visibility == "informed"
        frontier_count += frontier_hypothesis
        normalized_task: Dict[str, Any] = {
            "visibility": visibility,
            "focus": focus,
        }
        if methodology_lens:
            normalized_task["methodology_lens"] = methodology_lens
        if frontier_hypothesis:
            normalized_task["frontier_hypothesis"] = True
        normalized.append(normalized_task)

    if blind_count < 2 or informed_count < 3:
        raise ValueError("tasks require at least 2 blind and 3 informed entries")
    if len(normalized) > MIN_TASKS and not _normalized_text(expansion_reason):
        raise ValueError("expansion_reason required when task count exceeds 5")
    if frontier_count > 2:
        raise ValueError("tasks may contain at most 2 frontier hypotheses")
    if frontier_count == 2 and not _normalized_text(frontier_reason):
        raise ValueError("frontier_reason required when 2 frontier hypotheses are used")
    return normalized


def _response_contract(*, frontier_hypothesis: bool = False) -> Dict[str, Any]:
    contract = {
        "task_id": "copy the packet task_id exactly",
        "outcome": "direction | abstain",
        "analysis": "complete current-node chemical reasoning",
        "primary_direction": {
            "target_site": {
                "kind": "bond | functional_group",
                "atom_indices": [],
                "bond_index": None,
                "description": "current-molecule site",
            },
            "direction_summary": "strategic chemical change",
            "synthesis_logic": "forward synthesis logic",
            "precursor_logic": "conceptual direct-precursor logic",
            "realization_steps": [],
            "first_executable_step": "one real next event",
            "strategic_value": "why this direction matters",
            "risks": [],
            "uncertainties": [],
        },
        "abstain_reason": "required only when outcome=abstain",
    }
    if frontier_hypothesis:
        contract["analysis"] = (
            "host-retained frontier analysis satisfying frontier_analysis_contract"
        )
        contract["primary_direction"].update({
            "hypothesis_status": FRONTIER_HYPOTHESIS_STATUS,
            "hypothesis_basis": (
                "compact analogy, mechanistic basis, and speculative boundary"
            ),
        })
    return contract


def build_round(
    *,
    round_id: str,
    tasks: Any,
    expansion_reason: str,
    frontier_reason: str,
    session_id: str,
    tree_id: str,
    node_id: str,
    canonical_smiles: str,
    structure: Dict[str, Any],
    hard_constraints: Dict[str, Any],
    route_thesis: Dict[str, Any],
    guidance: List[Dict[str, Any]],
    parent_step: Dict[str, Any],
    current_action_space: Dict[str, Any],
    failure_summary: List[Dict[str, Any]],
    knowledge_profile_digest: str,
) -> Dict[str, Any]:
    """Build one immutable, unpersisted scouting round."""
    expansion_reason = _normalized_text(expansion_reason)
    frontier_reason = _normalized_text(frontier_reason)
    normalized_tasks = _normalize_tasks(
        tasks,
        expansion_reason=expansion_reason,
        frontier_reason=frontier_reason,
    )
    if not round_id or not session_id or not tree_id or not node_id or not canonical_smiles:
        raise ValueError("round binding requires session, tree, node, and structure identity")

    structure_projection = {
        "canonical_smiles": canonical_smiles,
        "atoms": copy.deepcopy(structure.get("atoms", [])),
        "bonds": copy.deepcopy(structure.get("bonds", [])),
        "rings": copy.deepcopy(structure.get("rings", [])),
        "stereo": copy.deepcopy(structure.get("stereo", {})),
    }
    response_contract = _response_contract()
    quality_contract = [
        "Return one primary direction: a bond disconnection or an FGI.",
        (
            "Name a concrete reaction idea and identify its current-structure target "
            "with packet indices; for example, give the named forward logic for "
            "disconnecting bond_index k at atoms [i, j], or name an FGI at its "
            "atom_indices."
        ),
        "Stop at reasonable direct precursors; do not generate a complete route.",
        "Keep chemistry uncertainty explicit and do not claim validation.",
    ]

    task_outputs: List[Dict[str, Any]] = []
    binding_tasks: List[Dict[str, Any]] = []
    for index, task in enumerate(normalized_tasks, 1):
        task_id = f"{round_id}:task:{index:02d}"
        task_response_contract = _response_contract(
            frontier_hypothesis=bool(task.get("frontier_hypothesis"))
        )
        packet: Dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "round_id": round_id,
            "task_id": task_id,
            "visibility": task["visibility"],
            "focus": task["focus"],
            "current_structure": copy.deepcopy(structure_projection),
            "hard_constraints": copy.deepcopy(hard_constraints),
            "quality_contract": list(quality_contract),
            "response_contract": task_response_contract,
        }
        if task.get("methodology_lens"):
            packet["methodology_lens"] = task["methodology_lens"]
            packet["methodology_lens_prompt"] = METHODOLOGY_LENSES[
                task["methodology_lens"]
            ]
        if task.get("frontier_hypothesis"):
            packet["frontier_hypothesis"] = True
            packet["frontier_analysis_contract"] = copy.deepcopy(
                FRONTIER_ANALYSIS_CONTRACT
            )
        if task["visibility"] == "informed":
            packet.update({
                "route_thesis": copy.deepcopy(route_thesis),
                "parent_step": copy.deepcopy(parent_step),
                "current_action_space": copy.deepcopy(current_action_space),
                "failure_summary": copy.deepcopy(failure_summary),
            })
        input_digest = _digest(packet)
        task_output = {
            "task_id": task_id,
            "visibility": task["visibility"],
            "focus": task["focus"],
            "input_digest": input_digest,
            "packet": packet,
        }
        binding_task = {
            "task_id": task_id,
            "visibility": task["visibility"],
            "focus": task["focus"],
            "input_digest": input_digest,
        }
        for key in ("methodology_lens", "frontier_hypothesis"):
            if key in task:
                task_output[key] = copy.deepcopy(task[key])
                binding_task[key] = copy.deepcopy(task[key])
        task_outputs.append(task_output)
        binding_tasks.append(binding_task)

    binding: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "round_id": round_id,
        "session_id": session_id,
        "tree_id": tree_id,
        "node_id": node_id,
        "canonical_smiles": canonical_smiles,
        "task_count": len(binding_tasks),
        "tasks": binding_tasks,
        "structure_digest": _digest(structure_projection),
        "hard_constraints_digest": _digest(hard_constraints),
        "route_plan_digest": _digest(route_thesis),
        "guidance_digest": _digest(guidance),
        "parent_step_digest": _digest(parent_step),
        "action_space_digest": _digest(current_action_space),
        "failure_summary_digest": _digest(failure_summary),
        "knowledge_profile_digest": knowledge_profile_digest,
    }
    if expansion_reason:
        binding["expansion_reason"] = expansion_reason
    if frontier_reason:
        binding["frontier_reason"] = frontier_reason
    binding["binding_digest"] = _digest(binding)

    return {
        "ok": True,
        "round_id": round_id,
        "round_binding": binding,
        "tasks": task_outputs,
        "response_contract": response_contract,
        "persisted": False,
        "next_step": "run scouts independently, then call scout_record",
    }

def _string_list(value: Any, field: str) -> List[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"{field} must be a list")
    normalized = [_normalized_text(item) for item in value]
    if any(not item for item in normalized):
        raise ValueError(f"{field} must not contain empty values")
    return normalized


def _normalize_target_site(
    raw: Any,
    structure: Dict[str, Any],
    *,
    field: str,
) -> Dict[str, Any]:
    if not isinstance(raw, dict):
        raise ValueError(f"{field} must be an object")
    kind = _normalized_text(raw.get("kind")).lower()
    if kind not in {"bond", "functional_group"}:
        raise ValueError(f"{field}.kind must be bond or functional_group")

    atom_indices = raw.get("atom_indices")
    if not isinstance(atom_indices, list) or not atom_indices:
        raise ValueError(f"{field}.atom_indices must be a non-empty list")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in atom_indices):
        raise ValueError(f"{field}.atom_indices must contain integers")
    if len(atom_indices) != len(set(atom_indices)):
        raise ValueError(f"{field}.atom_indices must be unique")

    valid_atom_indices = {
        atom.get("idx")
        for atom in structure.get("atoms", []) or []
        if isinstance(atom, dict)
    }
    if any(index not in valid_atom_indices for index in atom_indices):
        raise ValueError(f"{field}.atom_indices reference the wrong structure")

    bond_index = raw.get("bond_index")
    valid_bonds = {
        bond.get("idx"): bond
        for bond in structure.get("bonds", []) or []
        if isinstance(bond, dict)
    }
    if bond_index is not None and (
        isinstance(bond_index, bool)
        or not isinstance(bond_index, int)
        or bond_index not in valid_bonds
    ):
        raise ValueError(f"{field}.bond_index references the wrong structure")
    if kind == "bond":
        if bond_index is None:
            raise ValueError(f"{field}.bond_index required for bond direction")
        bond_atoms = set(valid_bonds[bond_index].get("atoms", []) or [])
        if set(atom_indices) != bond_atoms:
            raise ValueError(f"{field}.atom_indices do not match bond_index")
    elif bond_index is not None:
        raise ValueError(f"{field}.bond_index must be null for functional_group")

    description = _normalized_text(raw.get("description"))
    if not description:
        raise ValueError(f"{field}.description required")
    return {
        "kind": kind,
        "atom_indices": list(atom_indices),
        "bond_index": bond_index,
        "description": description,
    }


def _normalize_direction(
    raw: Any,
    structure: Dict[str, Any],
    *,
    field: str,
    frontier_hypothesis: bool = False,
) -> Dict[str, Any]:
    if not isinstance(raw, dict):
        raise ValueError(f"{field} must be an object")
    direction = {
        "target_site": _normalize_target_site(
            raw.get("target_site"),
            structure,
            field=f"{field}.target_site",
        )
    }
    for key in (
        "direction_summary",
        "synthesis_logic",
        "precursor_logic",
        "first_executable_step",
        "strategic_value",
    ):
        value = _normalized_text(raw.get(key))
        if not value:
            raise ValueError(f"{field}.{key} required")
        direction[key] = value
    direction["realization_steps"] = _string_list(
        raw.get("realization_steps"),
        f"{field}.realization_steps",
    )
    direction["risks"] = _string_list(raw.get("risks"), f"{field}.risks")
    direction["uncertainties"] = _string_list(
        raw.get("uncertainties"),
        f"{field}.uncertainties",
    )
    hypothesis_status = _normalized_text(raw.get("hypothesis_status")).lower()
    hypothesis_basis = _normalized_text(raw.get("hypothesis_basis"))
    if frontier_hypothesis:
        if hypothesis_status != FRONTIER_HYPOTHESIS_STATUS:
            raise ValueError(
                f"{field}.hypothesis_status must be {FRONTIER_HYPOTHESIS_STATUS}"
            )
        if not hypothesis_basis:
            raise ValueError(f"{field}.hypothesis_basis required")
    elif hypothesis_status == FRONTIER_HYPOTHESIS_STATUS:
        raise ValueError(
            f"{field}.hypothesis_status requires a frontier_hypothesis task"
        )
    elif hypothesis_status and hypothesis_status not in HYPOTHESIS_STATUSES:
        raise ValueError(f"{field}.hypothesis_status is invalid")
    if hypothesis_status:
        direction["hypothesis_status"] = hypothesis_status
    if hypothesis_basis:
        direction["hypothesis_basis"] = hypothesis_basis
    return direction


def review_projection(record: Dict[str, Any]) -> Dict[str, Any]:
    results = {
        item.get("task_id"): item
        for item in record.get("results", []) or []
        if isinstance(item, dict)
    }
    shortlist = []
    for task_id in record.get("shortlist_task_ids", []) or []:
        result = results.get(task_id) or {}
        if result.get("outcome") == "direction":
            shortlist.append({
                "task_id": task_id,
                "primary_direction": copy.deepcopy(result.get("primary_direction", {})),
            })
    return {
        "round_id": record.get("round_id", ""),
        "node_id": record.get("node_id", ""),
        "review_summary": record.get("review_summary", ""),
        "shortlist_task_ids": list(record.get("shortlist_task_ids", []) or []),
        "shortlisted_directions": shortlist,
    }


def build_prompt_view(
    rounds: List[Dict[str, Any]],
    node_id: str,
) -> Dict[str, Any]:
    """Project only the latest recorded review for the active node."""
    node_rounds = [
        record
        for record in rounds
        if isinstance(record, dict) and record.get("node_id") == node_id
    ]
    if not node_rounds:
        return {}
    return {
        "scouting_round_count": len(node_rounds),
        "scouting_review": review_projection(node_rounds[-1]),
        "scouting_guidance": SCOUTING_PROMPT_GUIDANCE,
    }


def build_report_view(
    rounds: List[Dict[str, Any]],
    tree: Any = None,
) -> Dict[str, Any]:
    """Build one read-only review and adoption projection for reports."""
    latest_records: Dict[str, Dict[str, Any]] = {}
    round_count_by_node: Dict[str, int] = {}
    seeds_by_node: Dict[str, List[Dict[str, Any]]] = {}
    directions_by_source: Dict[tuple[str, str], Dict[str, Any]] = {}
    for record in rounds:
        if not isinstance(record, dict):
            continue
        node_id = _normalized_text(record.get("node_id"))
        if not node_id:
            continue
        latest_records[node_id] = record
        round_count_by_node[node_id] = round_count_by_node.get(node_id, 0) + 1
        seeds = record.get("deferred_review_seeds", []) or []
        if seeds:
            seeds_by_node.setdefault(node_id, []).extend(copy.deepcopy(seeds))
        round_id = _normalized_text(record.get("round_id"))
        for result in record.get("results", []) or []:
            if not isinstance(result, dict):
                continue
            task_id = _normalized_text(result.get("task_id"))
            direction = result.get("primary_direction") or {}
            if round_id and task_id and isinstance(direction, dict):
                directions_by_source[(round_id, task_id)] = direction

    adopted_sources: List[Dict[str, Any]] = []
    for reaction in getattr(tree, "reaction_nodes", []) or []:
        decision = getattr(reaction, "llm_decision", None)
        audit = getattr(decision, "decision_audit", {}) or {}
        source = audit.get("scouting_source")
        if source:
            direction = directions_by_source.get((
                _normalized_text(source.get("round_id")),
                _normalized_text(source.get("task_id")),
            ), {})
            adopted_source = {
                "step_id": getattr(reaction, "step_id", ""),
                "node_id": getattr(reaction, "product_node", ""),
                "scouting_source": copy.deepcopy(source),
            }
            if direction.get("hypothesis_status") == FRONTIER_HYPOTHESIS_STATUS:
                adopted_source["hypothesis_status"] = FRONTIER_HYPOTHESIS_STATUS
                adopted_source["hypothesis_basis"] = direction.get(
                    "hypothesis_basis", ""
                )
            adopted_sources.append(adopted_source)

    if not latest_records and not adopted_sources:
        return {}
    return {
        "latest_review_by_node": {
            node_id: review_projection(record)
            for node_id, record in latest_records.items()
        },
        "round_count_by_node": round_count_by_node,
        "deferred_review_seeds": seeds_by_node,
        "adopted_sources": adopted_sources,
    }


def resolve_review_seed(
    rounds: List[Dict[str, Any]],
    seed_id: Any,
    node_id: str,
) -> Dict[str, Any] | None:
    """Resolve one self-contained seed snapshot before variant pruning."""
    normalized_seed_id = _normalized_text(seed_id)
    if not normalized_seed_id:
        return None
    for record in rounds:
        if not isinstance(record, dict):
            continue
        for seed in record.get("deferred_review_seeds", []) or []:
            if not isinstance(seed, dict) or seed.get("seed_id") != normalized_seed_id:
                continue
            if seed.get("node_id") != node_id:
                raise ValueError("review_seed_id must belong to the reviewed node")
            return copy.deepcopy(seed)
    raise ValueError("review_seed_id not found")


def prune_rounds_for_nodes(
    rounds: List[Dict[str, Any]],
    node_ids: set[str],
) -> Dict[str, Any]:
    """Remove complete rounds for reviewed or pruned nodes."""
    retained: List[Dict[str, Any]] = []
    removed_round_ids: List[str] = []
    for record in rounds:
        if isinstance(record, dict) and record.get("node_id") in node_ids:
            removed_round_ids.append(str(record.get("round_id", "")))
        else:
            retained.append(record)
    return {
        "rounds": retained,
        "removed_round_ids": [
            round_id for round_id in removed_round_ids if round_id
        ],
    }


def validate_and_build_record(
    *,
    provided_binding: Any,
    expected_round: Dict[str, Any],
    results: Any,
    review_summary: Any,
    shortlist_task_ids: Any,
    deferred_review_seeds: Any,
    existing_rounds: List[Dict[str, Any]],
    created_at: str,
) -> Dict[str, Any]:
    """Validate a complete batch and return a record without mutating state."""
    expected_binding = expected_round.get("round_binding") or {}
    if not isinstance(provided_binding, dict) or provided_binding != expected_binding:
        raise ValueError("round_binding does not match current node inputs")

    task_rows = expected_binding.get("tasks", []) or []
    expected_task_ids = [row.get("task_id") for row in task_rows]
    task_rows_by_id = {row.get("task_id"): row for row in task_rows}
    if not isinstance(results, list) or len(results) != len(expected_task_ids):
        raise ValueError("results must contain one outcome for every task")

    structure = expected_round["tasks"][0]["packet"]["current_structure"]
    normalized_by_task: Dict[str, Dict[str, Any]] = {}
    for index, raw in enumerate(results):
        if not isinstance(raw, dict):
            raise ValueError(f"results[{index}] must be an object")
        task_id = _normalized_text(raw.get("task_id"))
        if task_id not in expected_task_ids:
            raise ValueError(f"results[{index}].task_id is not part of this round")
        if task_id in normalized_by_task:
            raise ValueError("results contain duplicate task_id")
        outcome = _normalized_text(raw.get("outcome")).lower()
        normalized: Dict[str, Any] = {"task_id": task_id, "outcome": outcome}
        if outcome == "direction":
            normalized["primary_direction"] = _normalize_direction(
                raw.get("primary_direction"),
                structure,
                field=f"results[{index}].primary_direction",
                frontier_hypothesis=bool(
                    (task_rows_by_id.get(task_id) or {}).get("frontier_hypothesis")
                ),
            )
        elif outcome == "abstain":
            if raw.get("primary_direction") not in (None, {}):
                raise ValueError(
                    f"results[{index}].primary_direction is not allowed for abstain"
                )
            reason = _normalized_text(raw.get("abstain_reason"))
            if not reason:
                raise ValueError(f"results[{index}].abstain_reason required")
            normalized["abstain_reason"] = reason
        elif outcome == "unavailable":
            if raw.get("primary_direction") not in (None, {}):
                raise ValueError(
                    f"results[{index}].primary_direction is not allowed for unavailable"
                )
            reason = _normalized_text(raw.get("unavailable_reason"))
            if reason:
                normalized["unavailable_reason"] = reason
        else:
            raise ValueError(f"results[{index}].outcome is invalid")
        normalized_by_task[task_id] = normalized

    normalized_results = [normalized_by_task[task_id] for task_id in expected_task_ids]
    review_summary = _normalized_text(review_summary)
    if not review_summary:
        raise ValueError("review_summary required")

    if shortlist_task_ids is None:
        shortlist_task_ids = []
    if not isinstance(shortlist_task_ids, list) or len(shortlist_task_ids) > 3:
        raise ValueError("shortlist_task_ids must be a list with at most 3 entries")
    shortlist = [_normalized_text(task_id) for task_id in shortlist_task_ids]
    if len(shortlist) != len(set(shortlist)):
        raise ValueError("shortlist_task_ids must be unique")
    for task_id in shortlist:
        if (normalized_by_task.get(task_id) or {}).get("outcome") != "direction":
            raise ValueError("shortlist_task_ids may reference only direction outcomes")

    if deferred_review_seeds is None:
        deferred_review_seeds = []
    if not isinstance(deferred_review_seeds, list) or len(deferred_review_seeds) > 3:
        raise ValueError("deferred_review_seeds must be a list with at most 3 entries")
    normalized_seeds: List[Dict[str, Any]] = []
    round_id = expected_binding["round_id"]
    node_id = expected_binding["node_id"]
    for index, raw in enumerate(deferred_review_seeds, 1):
        if not isinstance(raw, dict):
            raise ValueError(f"deferred_review_seeds[{index - 1}] must be an object")
        sources = raw.get("source_task_ids")
        if not isinstance(sources, list) or not sources:
            raise ValueError(f"deferred_review_seeds[{index - 1}].source_task_ids required")
        source_ids = [_normalized_text(task_id) for task_id in sources]
        if len(source_ids) != len(set(source_ids)):
            raise ValueError("seed source_task_ids must be unique")
        if any(
            (normalized_by_task.get(task_id) or {}).get("outcome") != "direction"
            for task_id in source_ids
        ):
            raise ValueError("seed source_task_ids may reference only direction outcomes")
        seed: Dict[str, Any] = {
            "seed_id": f"{round_id}:seed:{index:02d}",
            "round_id": round_id,
            "node_id": node_id,
            "source_task_ids": source_ids,
            "target_site": _normalize_target_site(
                raw.get("target_site"),
                structure,
                field=f"deferred_review_seeds[{index - 1}].target_site",
            ),
        }
        for key in (
            "direction_summary",
            "synthesis_logic",
            "precursor_logic",
            "first_executable_step",
            "why_deferred",
            "revisit_trigger",
        ):
            value = _normalized_text(raw.get(key))
            if not value:
                raise ValueError(f"deferred_review_seeds[{index - 1}].{key} required")
            seed[key] = value
        frontier_directions = [
            (normalized_by_task[task_id].get("primary_direction") or {})
            for task_id in source_ids
            if (
                normalized_by_task[task_id].get("primary_direction") or {}
            ).get("hypothesis_status") == FRONTIER_HYPOTHESIS_STATUS
        ]
        if frontier_directions:
            seed["hypothesis_status"] = FRONTIER_HYPOTHESIS_STATUS
            bases = list(
                dict.fromkeys(
                    direction.get("hypothesis_basis", "")
                    for direction in frontier_directions
                    if direction.get("hypothesis_basis")
                )
            )
            seed["hypothesis_basis"] = " | ".join(bases)
        seed["constraints"] = _string_list(
            raw.get("constraints"),
            f"deferred_review_seeds[{index - 1}].constraints",
        )
        normalized_seeds.append(seed)

    record: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "round_id": round_id,
        "node_id": node_id,
        "canonical_smiles": expected_binding["canonical_smiles"],
        "round_binding": copy.deepcopy(expected_binding),
        "results": normalized_results,
        "review_summary": review_summary,
        "shortlist_task_ids": shortlist,
        "deferred_review_seeds": normalized_seeds,
    }
    record_digest = _digest(record)
    for existing in existing_rounds:
        if existing.get("round_id") != round_id:
            continue
        if existing.get("record_digest") == record_digest:
            return {"record": copy.deepcopy(existing), "idempotent": True}
        raise ValueError("round_id already recorded with different content")

    record["record_digest"] = record_digest
    record["created_at"] = created_at
    return {"record": record, "idempotent": False}

def validate_scouting_source(
    *,
    rounds: List[Dict[str, Any]],
    node_id: str,
    source: Any,
    action: Dict[str, Any],
) -> Dict[str, str] | None:
    """Validate attempt provenance before any sandbox mutation."""
    if source is None:
        return None
    if not isinstance(source, dict):
        raise ValueError("scouting_source must be an object")
    round_id = _normalized_text(source.get("round_id"))
    task_id = _normalized_text(source.get("task_id"))
    if not round_id or not task_id:
        raise ValueError("scouting_source requires round_id and task_id")

    record = next(
        (item for item in rounds if item.get("round_id") == round_id),
        None,
    )
    if record is None:
        raise ValueError("scouting_source round_id not found")
    if record.get("node_id") != node_id:
        raise ValueError("scouting_source must belong to the active node")
    task_result = next(
        (item for item in record.get("results", []) if item.get("task_id") == task_id),
        None,
    )
    if not task_result or task_result.get("outcome") != "direction":
        raise ValueError("scouting_source task_id must reference a direction")

    adoption_reason = _normalized_text(source.get("adoption_reason"))
    if task_id not in (record.get("shortlist_task_ids", []) or []) and not adoption_reason:
        raise ValueError("adoption_reason required for a non-shortlist scouting source")

    target_site = (task_result.get("primary_direction") or {}).get("target_site") or {}
    direction_kind = target_site.get("kind")
    action_source = action.get("source")
    if direction_kind == "bond" and action_source == "fgi":
        raise ValueError("scouting_source direction does not match the action kind")
    if direction_kind == "functional_group" and action_source == "bond":
        raise ValueError("scouting_source direction does not match the action kind")

    target_bond = target_site.get("bond_index")
    action_bond = action.get("actual_bond_idx")
    if action_bond is None:
        action_bond = action.get("bond_idx")
    # System FGI actions expose a stable option identity but no atom map.
    site_evidence_matched = (
        action_source == "fgi" and action.get("fgi_idx") is not None
    )
    if target_bond is not None and action_bond is not None:
        if int(target_bond) != int(action_bond):
            raise ValueError("scouting_source target bond does not match the action")
        site_evidence_matched = True
    target_atoms = set(target_site.get("atom_indices", []) or [])
    action_atoms = set(action.get("atoms", []) or [])
    if target_atoms and action_atoms and target_atoms != action_atoms:
        raise ValueError("scouting_source target atoms do not match the action")
    if target_atoms and action_atoms:
        site_evidence_matched = True
    changed_atom_sets = [
        set(change.get("product_atoms", []) or [])
        for change in action.get("changed_bonds", []) or []
        if isinstance(change, dict)
        and isinstance(change.get("product_atoms"), list)
        and change.get("product_atoms")
    ]
    if target_atoms and changed_atom_sets and target_atoms not in changed_atom_sets:
        raise ValueError("scouting_source target atoms do not match the action")
    if target_atoms and changed_atom_sets:
        site_evidence_matched = True
    if not site_evidence_matched:
        raise ValueError(
            "scouting_source action does not declare enough target-site evidence"
        )

    normalized = {"round_id": round_id, "task_id": task_id}
    if adoption_reason:
        normalized["adoption_reason"] = adoption_reason
    return normalized
