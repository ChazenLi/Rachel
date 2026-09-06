"""
BFS 逆合成编排引擎
==================
流程控制器，不做化学决策。所有判断信息传给 LLM，由 LLM 决定。

核心循环（分层上下文 + 沙盒试探）:
    prepare_next()       → compact 分子认知
    reaction_sites()     → 第一层 site/reaction 菜单
    explore_site(id)     → 第二层同 site 候选
    try_candidate(id)    → 标准化候选沙盒验证
    commit_decision()    → 确认满意，正式写入树
    accept_terminal()    → 标记 terminal
    skip_current()       → 跳过

安全阀:
    max_depth / max_steps / max_queue_size
"""

from __future__ import annotations

import json
import logging
import time
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

from Rachel.chem_tools._rdkit_utils import (
    canonical,
    has_assigned_stereocenter,
    parse_mol,
    validate_smiles,
)
from Rachel.chem_tools.cs_score import CS_TRIVIAL, compute_cs_score, classify_complexity
from Rachel.chem_tools.bond_break import (
    execute_disconnection, execute_fgi, preview_disconnections,
    try_retro_template,
)
from Rachel.chem_tools.forward_validate import validate_forward, check_atom_balance
from Rachel.chem_tools.fg_detect import detect_functional_groups
from Rachel.chem_tools.fg_warnings import suggest_protection_needs
from Rachel.chem_tools.mol_info import analyze_molecule
from Rachel.chem_tools.precursor_normalization import normalize_precursor_set
from Rachel.chem_tools.site_audit import audit_site_retention
from Rachel.chem_tools.topology_intent import action_declares_topology_change

from Rachel.tools.llm_retro_platform import (
    _compress_fg_instances,
    build_decision_context,
    format_context_text,
)

from .retro_tree import (
    RetrosynthesisTree,
    MoleculeNode,
    ReactionNode,
    MoleculeRole,
    TemplateEvidence,
    LLMDecision,
    parse_precursors,
)
from .retro_state import SynthesisAuditState
from .strategy_disclosure import (
    build_reaction_opportunity_brief,
    build_reaction_index,
    build_site_reaction_map,
    explore_reaction as expand_reaction_family,
    explore_site as expand_site_candidates,
)


def _collect_knowledge_refs(payload: Any) -> List[Dict[str, str]]:
    """Collect path-free knowledge references already present in a payload."""
    refs: List[Dict[str, str]] = []
    seen: Set[Tuple[Tuple[str, str], ...]] = set()

    def add_ref(value: Any) -> None:
        if not isinstance(value, dict):
            return
        ref = {
            str(key): str(item)
            for key, item in value.items()
            if isinstance(item, (str, int, float, bool))
        }
        if not ref or "pack_id" not in ref or "entry_id" not in ref:
            return
        key = tuple(sorted(ref.items()))
        if key not in seen:
            seen.add(key)
            refs.append(ref)

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                if key == "knowledge_ref":
                    add_ref(item)
                elif key == "knowledge_refs" and isinstance(item, list):
                    for ref in item:
                        add_ref(ref)
                else:
                    visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(payload)
    return refs
from .prompt_mount import build_prompt_brief, build_prompt_mount

logger = logging.getLogger(__name__)


def _build_validation_micro(forward_validation: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Keep only the gate-critical validation fields in sandbox/session payloads."""
    if not isinstance(forward_validation, dict):
        return {}

    checks = forward_validation.get("checks", {})
    assessment = forward_validation.get("assessment", {})
    evidence_packet = assessment.get("evidence_packet", {}) or {}
    score_breakdown = assessment.get("score_breakdown", {})
    template_exec = checks.get("template_execution", {})
    graph_delta = checks.get("graph_delta", {}) or {}
    atom_mapping = checks.get("atom_mapping", {}) or {}
    fg_delta = checks.get("fg_delta", {}) or {}
    ring_topology = checks.get("ring_topology", {}) or {}
    ring_family = checks.get("ring_family", {}) or {}
    family = ring_family.get("family_interpretation") or {}
    policy_preview = assessment.get("policy_preview") or {}
    gate = assessment.get("gate", {}) or {}
    policy_adjustments = gate.get("policy_adjustments", []) or []
    topology_codes = [
        str(item.get("code", ""))
        for item in (ring_topology.get("violations", []) or [])
        if item.get("code")
    ]
    topology_codes.extend(
        str(item.get("code", ""))
        for item in (ring_family.get("violations", []) or [])
        if item.get("code")
    )

    template_attempted = bool(template_exec.get("attempted", False))
    template_status = (
        "not_attempted"
        if not template_attempted
        else "target_regenerated"
        if template_exec.get("target_in_products") is True
        else "target_not_regenerated"
    )
    micro = {
        "template_attempted": template_attempted,
        "template_status": template_status,
        "scaffold_alignment": float(score_breakdown.get("scaffold_alignment", 0.0)),
        "graph_delta_risk": str(graph_delta.get("risk_level", "") or ""),
        "graph_delta_codes": list(graph_delta.get("delta_codes", []) or []),
        "atom_mapping_status": str(atom_mapping.get("status", "") or ""),
        "atom_mapping_confidence": str(atom_mapping.get("confidence", "") or ""),
        "atom_mapping_codes": list(atom_mapping.get("delta_codes", []) or []),
        "atom_mapping_formed_bond_count": len(atom_mapping.get("formed_bonds", []) or []),
        "atom_mapping_cleaved_bond_count": len(atom_mapping.get("cleaved_bonds", []) or []),
        "atom_mapping_mcs_product_coverage": atom_mapping.get("mcs_product_coverage"),
        "atom_mapping_ambiguous": bool((atom_mapping.get("ambiguity") or {}).get("ambiguous", False)),
        "atom_mapping_missing_fields": list(
            (atom_mapping.get("declared_evidence") or {}).get("missing_fields", []) or []
        ),
        "fg_delta_codes": list(fg_delta.get("delta_codes", []) or []),
        "ring_topology_risk": str(ring_topology.get("risk_level", "") or ""),
        "ring_delta_codes": list(ring_topology.get("delta_codes", []) or []),
        "topology_violation_codes": list(dict.fromkeys(topology_codes)),
        "family_key": str(family.get("family_key", "") or ""),
        "family_class": str(family.get("family_class", "") or ""),
        "family_state": str(family.get("state", "") or ""),
        "family_explained_delta_codes": list(family.get("explained_deltas", []) or []),
        "family_unexplained_delta_codes": list(family.get("unexplained_deltas", []) or []),
        "family_forbidden_delta_conflicts": list(family.get("forbidden_delta_conflicts", []) or []),
        "family_supporting_codes": list(family.get("supporting_codes", []) or []),
        "family_required_evidence": list(family.get("required_evidence", []) or []),
        "family_evidence_state": str((family.get("family_evidence") or {}).get("state", "") or ""),
        "family_provided_evidence_keys": list(
            (family.get("family_evidence") or {}).get("provided_keys", []) or []
        ),
        "family_missing_evidence_keys": list(
            (family.get("family_evidence") or {}).get("missing_keys", []) or []
        ),
        "family_missing_evidence_codes": list(
            (family.get("family_evidence") or {}).get("missing_evidence_codes", []) or []
        ),
        "policy_preview": {
            key: value
            for key, value in {
                "dry_run_only": policy_preview.get("dry_run_only"),
                "policy_source": policy_preview.get("policy_source", ""),
                "family_key": policy_preview.get("family_key", ""),
                "family_class": policy_preview.get("family_class", ""),
                "family_state": policy_preview.get("family_state", ""),
                "risk_level": policy_preview.get("risk_level", ""),
                "current_gate_state": policy_preview.get("current_gate_state", ""),
                "current_commit_policy": policy_preview.get("current_commit_policy", ""),
                "preview_gate_state": policy_preview.get("preview_gate_state", ""),
                "preview_commit_policy": policy_preview.get("preview_commit_policy", ""),
                "applied_downgrade_codes": list(policy_preview.get("applied_downgrade_codes", []) or []),
                "preserved_gate_codes": list(policy_preview.get("preserved_gate_codes", []) or []),
                "preview_hard_block_codes": list(policy_preview.get("preview_hard_block_codes", []) or []),
                "preview_override_allowed_codes": list(policy_preview.get("preview_override_allowed_codes", []) or []),
                "preview_missing_evidence_codes": list(policy_preview.get("preview_missing_evidence_codes", []) or []),
                "preview_soft_warning_codes": list(policy_preview.get("preview_soft_warning_codes", []) or []),
                "would_change_gate": policy_preview.get("would_change_gate"),
                "would_change_commit_policy": policy_preview.get("would_change_commit_policy"),
            }.items()
            if value not in (None, "", [], {})
        },
        "policy_adjustments": [
            {
                key: value
                for key, value in {
                    "activation_scope": item.get("activation_scope", ""),
                    "policy_source": item.get("policy_source", ""),
                    "family_key": item.get("family_key", ""),
                    "family_class": item.get("family_class", ""),
                    "family_state": item.get("family_state", ""),
                    "risk_level": item.get("risk_level", ""),
                    "from_bucket": item.get("from_bucket", ""),
                    "to_bucket": item.get("to_bucket", ""),
                    "applied_downgrade_codes": list(item.get("applied_downgrade_codes", []) or []),
                }.items()
                if value not in (None, "", [], {})
            }
            for item in policy_adjustments
            if isinstance(item, dict)
        ],
        "evidence_packet": evidence_packet,
    }
    if template_attempted:
        micro["template_product_similarity"] = float(
            score_breakdown.get("template_match", 0.0)
        )
    return micro


def _validation_unavailable_payload(
    exc: Exception,
    precursors: List[str],
    product_smiles: str,
    action_context: Dict[str, Any],
) -> Dict[str, Any]:
    """Represent validator failures as an explicit non-committable gate."""
    exception_type = exc.__class__.__name__
    message = str(exc) or exception_type
    validation_error = {
        "code": "validation_unavailable",
        "exception_type": exception_type,
        "message": message,
    }
    gate = {
        "gate_state": "hard_block",
        "commit_policy": "block",
        "hard_blocks": [
            {
                "code": "validation_unavailable",
                "category": "validation",
                "message": "forward validation failed before producing a gate",
                "evidence": validation_error,
            }
        ],
        "soft_warnings": [],
        "override_allowed": [],
        "missing_evidence": [],
        "llm_override_allowed": False,
        "recommended_action": "Rerun validation or choose another action before commit.",
    }
    policy_preview = {
        "policy_source": "validation_unavailable",
        "current_gate_state": "hard_block",
        "current_commit_policy": "block",
        "preview_gate_state": "hard_block",
        "preview_commit_policy": "block",
        "preview_hard_block_codes": ["validation_unavailable"],
        "preview_override_allowed_codes": [],
        "preview_missing_evidence_codes": [],
        "preview_soft_warning_codes": [],
        "would_change_gate": False,
        "would_change_commit_policy": False,
    }
    evidence_packet = {
        "validation_error": validation_error,
        "action_context": dict(action_context or {}),
        "precursors": list(precursors or []),
        "product_smiles": product_smiles,
        "recommended_action": gate["recommended_action"],
    }
    return {
        "ok": False,
        "precursors": list(precursors or []),
        "product_smiles": product_smiles,
        "checks": {
            "action_context": dict(action_context or {}),
            "validation_error": validation_error,
        },
        "assessment": {
            "pass": False,
            "feasibility_score": 0.0,
            "hard_fail_reasons": ["validation_unavailable"],
            "gate": gate,
            "policy_preview": policy_preview,
            "evidence_packet": evidence_packet,
        },
    }


def _custom_context_implies_ring_action(action_context: Dict[str, Any], reaction_type: str) -> bool:
    """Return explicit custom topology intent without reaction-name guessing."""
    del reaction_type
    return action_declares_topology_change(action_context)


def _invalid_precursor_smiles(precursors: List[str]) -> List[Dict[str, str]]:
    """Return invalid precursor SMILES with strict RDKit/sanitize reasons."""
    invalid: List[Dict[str, str]] = []
    for smi in precursors:
        ok, reason = validate_smiles(smi)
        if not ok:
            invalid.append({"smiles": smi, "reason": reason})
    return invalid


def _preflight_precursors(precursors: List[str]) -> Dict[str, Any]:
    """Canonicalize precursors and record organometallic source obligations."""
    result = normalize_precursor_set(precursors)
    return {
        "precursors": list(result.get("precursors", []) or []),
        "normalizations": list(result.get("normalizations", []) or []),
        "invalid": list(result.get("invalid", []) or []),
        "changed": bool(result.get("changed", False)),
    }


def _invalid_precursor_detail(invalid_precursors: List[Dict[str, Any]]) -> str:
    return "; ".join(
        f"{item.get('smiles', '')} ({item.get('reason', '')})"
        for item in invalid_precursors
    )


def _precursor_normalization_payload(
    original_precursors: List[str],
    preflight: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    normalizations = list(preflight.get("normalizations", []) or [])
    if not normalizations:
        return None
    return {
        "original_precursors": list(original_precursors),
        "normalized_precursors": list(preflight.get("precursors", []) or []),
        "normalizations": normalizations,
        "invalid": list(preflight.get("invalid", []) or []),
        "changed": bool(preflight.get("changed", False)),
    }


def _build_bond_fg_context(smiles: str, bond_atoms: Tuple[int, int]) -> Dict[str, List[str]]:
    """Return a local FG micro-audit for explore_bond without expanding compact context."""
    empty_context = {"atom_i": [], "atom_j": [], "spanning": [], "nearby": []}
    if len(bond_atoms) != 2:
        return empty_context

    mol = parse_mol(smiles)
    fg_result = detect_functional_groups(smiles)
    groups_data = fg_result.get("groups", {})
    if mol is None or not isinstance(groups_data, dict):
        return empty_context

    try:
        atom_i, atom_j = int(bond_atoms[0]), int(bond_atoms[1])
        local_shell = {atom_i, atom_j}
        for atom_idx in (atom_i, atom_j):
            atom = mol.GetAtomWithIdx(atom_idx)
            local_shell.update(nb.GetIdx() for nb in atom.GetNeighbors())
    except Exception:
        return empty_context

    # Recompute the local FG view on demand from the raw detector output. This keeps
    # the complete SMARTS hit map in the tool layer and exposes only a bond-specific
    # micro-audit at explore time, which is the right layer for site verification.
    compressed_instances = _compress_fg_instances(groups_data)
    seen: Dict[str, Set[str]] = {key: set() for key in empty_context}
    bond_fg_context: Dict[str, List[str]] = {key: [] for key in empty_context}
    for inst in compressed_instances:
        display_name = str(inst.get("display_name") or inst.get("name") or "")
        atom_set = set(inst.get("atoms", []))
        if not display_name or not atom_set:
            continue

        if atom_i in atom_set and atom_j in atom_set:
            bucket = "spanning"
        elif atom_i in atom_set:
            bucket = "atom_i"
        elif atom_j in atom_set:
            bucket = "atom_j"
        elif atom_set & local_shell:
            bucket = "nearby"
        else:
            continue

        if display_name not in seen[bucket]:
            seen[bucket].add(display_name)
            bond_fg_context[bucket].append(display_name)

    return bond_fg_context


# ─────────────────────────────────────────────────────────────────────────
# 编排器数据类
# ─────────────────────────────────────────────────────────────────────────

@dataclass
class ProposalContext:
    """prepare_next 的返回值：当前分子的分析 + 断键提案，供 LLM 决策。

    分层上下文策略:
      - to_dict()  → 精简摘要：分子信息 + top-N 键位概览（不含完整前体列表）
      - to_dict(detail="full") → 完整版（兼容旧行为）
      - explore_bond(idx) / explore_fgi() → 按需展开细节（由编排器代理调用）
    """
    smiles: str
    node_id: str
    depth: int
    cs_score: float = 0.0
    classification: str = ""
    is_terminal: bool = False
    is_target: bool = False
    depth_limited: bool = False
    decision_context: Optional[Dict[str, Any]] = None
    seen_smiles: Set[str] = field(default_factory=set)
    steps_executed: int = 0
    steps_remaining: int = 0
    queue_preview: List[Dict[str, Any]] = field(default_factory=list)
    audit_state_summary: Dict[str, Any] = field(default_factory=dict)
    failed_attempts_for_current: List[Dict[str, Any]] = field(default_factory=list)
    decision_tier: str = "standard"     # quick_pass / standard
    knowledge_profile: Any = field(default=None, repr=False, compare=False)

    def _normalize_compact_window(
        self,
        *,
        top_n: int = 5,
        bond_offset: int = 0,
        bond_limit: Optional[int] = None,
        fgi_offset: int = 0,
        fgi_limit: int = 5,
    ) -> Tuple[int, int, int, int]:
        """Normalize compact window parameters while preserving the old top_n contract."""
        # Historical contract kept instead of being removed:
        # top_n used to control the default compact truncation window by itself.
        # We keep it alive for backward compatibility, but the new explicit window
        # parameters now define which compact slice is being shown.
        effective_bond_limit = top_n if bond_limit is None else bond_limit

        try:
            bond_offset = max(0, int(bond_offset))
        except (TypeError, ValueError):
            bond_offset = 0
        try:
            effective_bond_limit = max(0, int(effective_bond_limit))
        except (TypeError, ValueError):
            effective_bond_limit = max(0, int(top_n))
        try:
            fgi_offset = max(0, int(fgi_offset))
        except (TypeError, ValueError):
            fgi_offset = 0
        try:
            fgi_limit = max(0, int(fgi_limit))
        except (TypeError, ValueError):
            fgi_limit = 5

        return bond_offset, effective_bond_limit, fgi_offset, fgi_limit

    def _build_bond_fgi_summary_payload(
        self,
        *,
        top_n: int = 5,
        bond_offset: int = 0,
        bond_limit: Optional[int] = None,
        fgi_offset: int = 0,
        fgi_limit: int = 5,
    ) -> Dict[str, Any]:
        """Build legacy/debug bond and FGI windows for diagnostic context only."""
        payload: Dict[str, Any] = {}
        if not self.decision_context:
            return payload

        bond_offset, bond_limit, fgi_offset, fgi_limit = self._normalize_compact_window(
            top_n=top_n,
            bond_offset=bond_offset,
            bond_limit=bond_limit,
            fgi_offset=fgi_offset,
            fgi_limit=fgi_limit,
        )

        ctx = self.decision_context
        bonds = ctx.get("disconnectable_bonds", [])
        bond_summary = []
        bond_end = min(len(bonds), bond_offset + bond_limit)
        for global_bond_idx in range(bond_offset, bond_end):
            b = bonds[global_bond_idx]
            alts = b.get("alternatives", [])
            summary = {
                "bond_idx": global_bond_idx,
                "actual_bond_idx": b.get("actual_bond_idx", -1),
                "atoms": b.get("atoms", []),
                "role_pair": b.get("role_pair", ["unclear", "unclear"]),
                "bond_type": b.get("bond_type", ""),
                "in_ring": b.get("in_ring", False),
                "heuristic_score": b.get("heuristic_score", 0),
                "n_alternatives": len(alts),
                "reaction_types": [
                    a.get("template", "").split("(")[0].strip()
                    for a in alts[:3]
                ],
            }
            sc = b.get("smart_capping", [])
            if sc:
                summary["smart_capping"] = [
                    {"type": c["reaction_type"], "conf": c["confidence"]}
                    for c in sc[:2]
                ]
            bond_summary.append(summary)

        payload["bond_summary"] = bond_summary
        payload["bond_summary_offset"] = bond_offset
        payload["bond_summary_limit"] = bond_limit
        if bond_end < len(bonds):
            payload["bonds_omitted"] = len(bonds) - bond_end

        fgi_list = ctx.get("fgi_options", [])
        if fgi_list:
            fgi_end = min(len(fgi_list), fgi_offset + fgi_limit)
            payload["fgi_summary"] = [
                {"fgi_idx": global_fgi_idx, "template": fgi_list[global_fgi_idx].get("template", "")}
                for global_fgi_idx in range(fgi_offset, fgi_end)
            ]
            payload["fgi_summary_offset"] = fgi_offset
            payload["fgi_summary_limit"] = fgi_limit
            if fgi_end < len(fgi_list):
                payload["fgi_omitted"] = len(fgi_list) - fgi_end

        return payload

    def _build_compact_payload(
        self,
        *,
        top_n: int = 5,
        bond_offset: int = 0,
        bond_limit: Optional[int] = None,
        fgi_offset: int = 0,
        fgi_limit: int = 5,
    ) -> Dict[str, Any]:
        """Build the default LLM cognition payload without the full first layer."""
        compact: Dict[str, Any] = {}
        if not self.decision_context:
            return compact

        ctx = self.decision_context
        molecule = ctx.get("molecule", {}) or {}
        compact["molecule_brief"] = {
            key: molecule[key]
            for key in (
                "formula",
                "mw",
                "heavy_atoms",
                "rings",
                "aromatic_rings",
                "stereocenters",
                "rotatable_bonds",
            )
            if key in molecule
        }
        compact["functional_group_brief"] = [
            {"name": fg.get("name", ""), "count": fg.get("count", 1)}
            for fg in ctx.get("functional_groups", []) or []
            if fg.get("name", "")
        ]
        complexity = ctx.get("complexity", {}) or {}
        compact["complexity_brief"] = {
            key: complexity[key]
            for key in ("cs_score", "classification", "is_terminal")
            if key in complexity
        }
        compact["n_bonds"] = ctx.get("n_bonds", 0)
        compact["n_fgi"] = len(ctx.get("fgi_options", []))
        compact["warnings"] = ctx.get("warnings", [])
        compact["primary_disclosure"] = "reaction_opportunity_brief"
        compact["reaction_opportunity_brief"] = build_reaction_opportunity_brief(ctx)
        compact["commands"] = {
            "first_layer": "reaction_sites()",
            "second_layer": "explore_site(site_id)",
            "structure_detail": 'context(detail="structure")',
            "custom_structure_detail": 'inspect_structures(molecules=[...], include_current=true)',
            "sandbox": "try_action(action_id)",
            "custom_action": "propose_action(...)",
            "terminal": "accept_terminal(reason=...)",
        }
        compact["hint"] = (
            "Use reaction_opportunity_brief for molecule-level triage; call "
            "reaction_sites() for the full first-layer site menu, then "
            "explore_site(site_id) and try_action(action_id)."
        )
        compact["prompt_brief"] = build_prompt_brief(
            build_prompt_mount(
                "context_compact",
                decision_context=ctx,
                knowledge_profile=self.knowledge_profile,
            )
        )
        return compact
    def _build_diagnostic_payload(
        self,
        *,
        top_n: int = 5,
        bond_offset: int = 0,
        bond_limit: Optional[int] = None,
        fgi_offset: int = 0,
        fgi_limit: int = 5,
    ) -> Dict[str, Any]:
        """Build debug/legacy context. This is not the default LLM payload."""
        diagnostic = self._build_compact_payload(
            top_n=top_n,
            bond_offset=bond_offset,
            bond_limit=bond_limit,
            fgi_offset=fgi_offset,
            fgi_limit=fgi_limit,
        )
        if not self.decision_context:
            return diagnostic

        reaction_index = build_reaction_index(self.decision_context)
        diagnostic["diagnostic_view"] = True
        diagnostic["reaction_index"] = reaction_index
        diagnostic["reaction_index_count"] = len(reaction_index)
        diagnostic.update(
            self._build_bond_fgi_summary_payload(
                top_n=top_n,
                bond_offset=bond_offset,
                bond_limit=bond_limit,
                fgi_offset=fgi_offset,
                fgi_limit=fgi_limit,
            )
        )
        return diagnostic

    def to_dict(
        self,
        detail: str = "compact",
        top_n: int = 5,
        *,
        bond_offset: int = 0,
        bond_limit: Optional[int] = None,
        fgi_offset: int = 0,
        fgi_limit: int = 5,
    ) -> Dict[str, Any]:
        """分层输出。

        detail:
          "compact" — 默认，键位只给概览（atoms, score, reaction_types, n_alternatives）
          "full"    — 完整版，包含所有前体方案（兼容旧行为，慎用）
        top_n: compact 模式下只展示前 N 个键位
        """
        if self.decision_tier == "quick_pass":
            quick = {
                "action": "awaiting_decision",
                "decision_tier": "quick_pass",
                "smiles": self.smiles,
                "node_id": self.node_id,
                "depth": self.depth,
                "cs_score": self.cs_score,
                "classification": self.classification,
                "is_terminal": self.is_terminal,
                "steps_executed": self.steps_executed,
                "hint": "Fast path: call accept to mark this molecule as terminal.",
            }
            if self.knowledge_profile is not None:
                quick["knowledge_profile_hash"] = self.knowledge_profile.digest
                quick["knowledge_refs"] = []
            return quick

        d: Dict[str, Any] = {
            "action": "awaiting_decision",
            "decision_tier": "standard",
            "smiles": self.smiles,
            "node_id": self.node_id,
            "depth": self.depth,
            "cs_score": self.cs_score,
            "classification": self.classification,
            "is_terminal": self.is_terminal,
            "is_target": self.is_target,
            "depth_limited": self.depth_limited,
            "steps_executed": self.steps_executed,
            "steps_remaining": self.steps_remaining,
        }

        if self.decision_context:
            ctx = self.decision_context
            if detail == "structure":
                d["molecule_structure"] = analyze_molecule(self.smiles)
            elif detail == "full":
                d["molecule"] = ctx.get("molecule", {})
                d["functional_groups"] = ctx.get("functional_groups", [])
                d["complexity"] = ctx.get("complexity", {})
                d["n_bonds"] = ctx.get("n_bonds", 0)
                d["n_fgi"] = len(ctx.get("fgi_options", []))
                d["warnings"] = ctx.get("warnings", [])
                d["disconnectable_bonds"] = ctx.get("disconnectable_bonds", [])
                d["fgi_options"] = ctx.get("fgi_options", [])
            elif detail == "diagnostic":
                d.update(
                    self._build_diagnostic_payload(
                        top_n=top_n,
                        bond_offset=bond_offset,
                        bond_limit=bond_limit,
                        fgi_offset=fgi_offset,
                        fgi_limit=fgi_limit,
                    )
                )
            else:
                # Historical inline compact builder kept as comments instead of being
                # deleted. Reason: compact is now explicitly windowed, and the session
                # path also reuses this exact builder so live/context/session stay in
                # one contract instead of drifting.
                # d["molecule"] = ctx.get("molecule", {})
                # d["functional_groups"] = ctx.get("functional_groups", [])
                # d["complexity"] = ctx.get("complexity", {})
                # d["n_bonds"] = ctx.get("n_bonds", 0)
                # d["n_fgi"] = len(ctx.get("fgi_options", []))
                # d["warnings"] = ctx.get("warnings", [])
                # bonds = ctx.get("disconnectable_bonds", [])
                # bond_summary = []
                # for i, b in enumerate(bonds[:top_n]):
                #     ...
                d.update(
                    self._build_compact_payload(
                        top_n=top_n,
                        bond_offset=bond_offset,
                        bond_limit=bond_limit,
                        fgi_offset=fgi_offset,
                        fgi_limit=fgi_limit,
                    )
                )

        if self.queue_preview:
            d["queue_preview"] = self.queue_preview
        if self.audit_state_summary:
            d["audit_state_summary"] = self.audit_state_summary
        if self.failed_attempts_for_current:
            d["failed_attempts_for_current"] = self.failed_attempts_for_current
        if self.knowledge_profile is not None:
            d["knowledge_profile_hash"] = self.knowledge_profile.digest
            d["knowledge_refs"] = _collect_knowledge_refs(
                [self.decision_context or {}, d]
            )
        return d

    def to_text(
        self,
        target_smiles: str = "",
        detail: str = "compact",
        *,
        bond_offset: int = 0,
        bond_limit: Optional[int] = None,
        fgi_offset: int = 0,
        fgi_limit: int = 5,
    ) -> str:
        """格式化为 LLM 可读文本。"""
        if self.decision_tier == "quick_pass":
            return (
                f"Molecule {self.smiles} - CS={self.cs_score:.2f} "
                f"({self.classification}); suitable as a terminal starting material."
            )
        parts = []
        parts.append(f"=== Retrosynthesis Decision - Step {self.steps_executed + 1} ===")
        if target_smiles:
            parts.append(f"Target product: {target_smiles}")
        parts.append(f"Current molecule: {self.smiles}  depth={self.depth}")
        parts.append(f"Complexity: CS={self.cs_score:.2f} ({self.classification})")
        parts.append("")

        if detail == "full" and self.decision_context:
            parts.append(format_context_text(self.decision_context))
        elif self.decision_context:
            ctx_payload = (
                self._build_diagnostic_payload(
                    top_n=5,
                    bond_offset=bond_offset,
                    bond_limit=bond_limit,
                    fgi_offset=fgi_offset,
                    fgi_limit=fgi_limit,
                )
                if detail == "diagnostic"
                else self._build_compact_payload(
                    top_n=5,
                    bond_offset=bond_offset,
                    bond_limit=bond_limit,
                    fgi_offset=fgi_offset,
                    fgi_limit=fgi_limit,
                )
            )
            fgs = ctx_payload.get("functional_group_brief", [])
            if fgs:
                fg_names = [g["name"] for g in fgs[:10]]
                parts.append(f"Functional groups: {', '.join(fg_names)}")

            opportunity = ctx_payload.get("reaction_opportunity_brief", {})
            if opportunity:
                parts.append(
                    "Reaction opportunity brief: "
                    f"{opportunity.get('site_count', 0)} sites, "
                    f"{opportunity.get('total_reaction_count', 0)} reactions, "
                    f"{opportunity.get('competing_site_count', 0)} competing sites"
                )
                high_sites = opportunity.get("high_competition_sites", []) or []
                for site in high_sites:
                    names = ", ".join(site.get("reaction_names", []) or [])
                    parts.append(
                        f"  {site.get('site_id', '')} | {site.get('site_hint', '')} "
                        f"({site.get('reaction_count', 0)} reactions: {names})"
                    )

            if detail == "diagnostic":
                bond_summary = ctx_payload.get("bond_summary", [])
                if bond_summary:
                    parts.append(f"Diagnostic bond window: {len(bond_summary)} entries")
                    for b in bond_summary:
                        role_pair = "/".join(b.get("role_pair", ["unclear", "unclear"]))
                        parts.append(
                            f"  [{b.get('bond_idx', -1)}] atoms={b.get('atoms', [])} "
                            f"rdkit_bond={b.get('actual_bond_idx', -1)} roles={role_pair}"
                        )
                fgi_summary = ctx_payload.get("fgi_summary", [])
                if fgi_summary:
                    parts.append(f"Diagnostic FGI window: {len(fgi_summary)} entries")
                    for f in fgi_summary:
                        parts.append(f"  [{f.get('fgi_idx', -1)}] {f.get('template', '')}")

            parts.append(
                "Hint: use reaction_sites() -> explore_site(site_id) -> "
                "try_action(action_id)."
            )
        return "\n".join(parts)


@dataclass
class CommitResult:
    """commit_decision 的返回值。"""
    success: bool
    reaction_node: Optional[ReactionNode] = None
    new_pending: List[str] = field(default_factory=list)
    new_terminal: List[str] = field(default_factory=list)
    tree_complete: bool = False
    cycle_warnings: List[str] = field(default_factory=list)
    forward_validation: Optional[Dict[str, Any]] = None
    knowledge_profile_hash: str = ""
    knowledge_refs: List[Dict[str, str]] = field(default_factory=list)
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {"success": self.success}
        if self.reaction_node:
            d["step_id"] = self.reaction_node.step_id
            d["reaction_smiles"] = self.reaction_node.reaction_smiles
        d["new_pending"] = self.new_pending
        d["new_terminal"] = self.new_terminal
        d["tree_complete"] = self.tree_complete
        if self.cycle_warnings:
            d["cycle_warnings"] = self.cycle_warnings
        if self.forward_validation:
            from .retro_tree import _flatten_fv
            d["forward_validation"] = _flatten_fv(self.forward_validation)
        if self.knowledge_profile_hash:
            d["knowledge_profile_hash"] = self.knowledge_profile_hash
        if self.knowledge_refs:
            d["knowledge_refs"] = self.knowledge_refs
        if self.error:
            d["error"] = self.error
        return d


@dataclass
class SandboxResult:
    """沙盒试探的返回值。不写入树，只返回前体 + 验证结果。

    LLM 可以反复调用 try_* 方法，比较不同方案，满意后再 commit。
    """
    success: bool
    precursors: List[str] = field(default_factory=list)
    reagents: List[str] = field(default_factory=list)
    precursor_details: List[Dict[str, Any]] = field(default_factory=list)
    forward_validation: Optional[Dict[str, Any]] = None
    atom_balance: Optional[Dict[str, Any]] = None
    cycle_warnings: List[str] = field(default_factory=list)
    reaction_type: str = ""
    template_id: str = ""
    template_name: str = ""
    execution_mode: str = ""
    declared_template_id: str = ""
    executed_template_id: str = ""
    candidate_consistency: Optional[Dict[str, Any]] = None
    precursor_normalization: Optional[Dict[str, Any]] = None
    validation_micro: Optional[Dict[str, Any]] = None
    site_audit: Optional[Dict[str, Any]] = None
    eas_site_audit: Optional[Dict[str, Any]] = None
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {"success": self.success}
        if self.precursors:
            d["precursors"] = self.precursors
        if self.reagents:
            d["reagents"] = self.reagents
        if self.precursor_details:
            d["precursor_details"] = self.precursor_details
        if self.forward_validation:
            from .retro_tree import _flatten_fv
            d["forward_validation"] = _flatten_fv(self.forward_validation)
        if self.validation_micro:
            d["validation_micro"] = self.validation_micro
        if self.atom_balance:
            d["atom_balance"] = {
                "balanced": self.atom_balance.get("balanced", False),
                "message": self.atom_balance.get("message", ""),
            }
        if self.cycle_warnings:
            d["cycle_warnings"] = self.cycle_warnings
        if self.reaction_type:
            d["reaction_type"] = self.reaction_type
        if self.template_id:
            d["template_id"] = self.template_id
        if self.template_name:
            d["template_name"] = self.template_name
        if self.execution_mode:
            d["execution_mode"] = self.execution_mode
        if self.declared_template_id:
            d["declared_template_id"] = self.declared_template_id
        if self.executed_template_id:
            d["executed_template_id"] = self.executed_template_id
        if self.candidate_consistency:
            d["candidate_consistency"] = self.candidate_consistency
        if self.precursor_normalization:
            d["precursor_normalization"] = self.precursor_normalization
        if self.site_audit:
            d["site_audit"] = self.site_audit
        if self.eas_site_audit:
            d["eas_site_audit"] = self.eas_site_audit
        if self.error:
            d["error"] = self.error
        d["hint"] = (
            "This sandbox result has not been written to the tree. "
            "If it is satisfactory, call commit."
        )
        return d


# ─────────────────────────────────────────────────────────────────────────
# 编排器
# ─────────────────────────────────────────────────────────────────────────

class RetrosynthesisOrchestrator:
    """BFS 驱动的逆合成编排引擎。

    编排器不做化学决策。is_terminal、depth_limited 等判断信息
    全部传给 LLM，由 LLM 决定是否继续拆解。

    终止判定三层:
      1. 重原子 ≤ 6 → 无条件 terminal
      2. CS score ≤ terminal_threshold → terminal
      3. 无可断键位 → terminal

    LLM 可覆盖: 对任何分子调用 accept_terminal()。
    """

    def __init__(
        self,
        target_smiles: str,
        target_name: str = "",
        *,
        max_depth: int = 15,
        max_steps: int = 50,
        max_queue_size: int = 200,
        terminal_cs_threshold: float = CS_TRIVIAL,
        auto_forward_validate: bool = True,
        knowledge_profile=None,
    ):
        can = canonical(target_smiles)
        if not can:
            raise ValueError(f"Invalid target SMILES: {target_smiles}")

        self.tree = RetrosynthesisTree(can, target_name)
        self.max_depth = max_depth
        self.max_steps = max_steps
        self.max_queue_size = max_queue_size
        self.terminal_cs_threshold = terminal_cs_threshold
        self.auto_forward_validate = auto_forward_validate
        if knowledge_profile is None:
            from Rachel.knowledge import get_base_profile

            knowledge_profile = get_base_profile()
        self.knowledge_profile = knowledge_profile

        # BFS 队列: (smiles, depth)
        self._queue: deque = deque()
        self._queue.append((can, 0))

        # 已见分子集合
        self._seen: Set[str] = {can}

        # 当前活跃 context
        self._current_context: Optional[ProposalContext] = None
        self._force_standard_smiles: Set[str] = set()

        # 统计
        self._steps_executed: int = 0
        self._start_time: float = time.time()

        # 审计状态
        self.audit_state = SynthesisAuditState()

        # 分析目标分子复杂度
        cs = compute_cs_score(can, knowledge_profile=self.knowledge_profile)
        self.tree.update_complexity(can, cs)
        self.audit_state.set_target_complexity(cs.get("cs_score", 0))

    # ── 状态查询 ──

    @property
    def pending_count(self) -> int:
        return len(self._queue)

    def step_budget_exhausted(self) -> bool:
        """Return whether unresolved work remains after the step budget is spent."""
        return self._steps_executed >= self.max_steps and bool(
            self._current_context is not None or self._queue
        )

    def _pending_precursors_for_queue(self, precursors: List[str]) -> List[str]:
        """Return unique nonterminal precursors that need new queue entries."""
        expanded_products = {
            reaction.product_node for reaction in self.tree.reaction_nodes
        }
        queued_smiles = {queued for queued, _ in self._queue}
        current_smiles = (
            self._current_context.smiles
            if self._current_context is not None
            else ""
        )
        pending: List[str] = []
        pending_set: Set[str] = set()
        for smiles in precursors:
            can = canonical(smiles) or smiles
            if can in pending_set:
                continue
            is_terminal, _ = self._check_terminal(can)
            if is_terminal:
                continue
            node = self.tree.get_molecule_by_smiles(can)
            already_expanded = bool(
                node and node.node_id in expanded_products
            )
            already_scheduled = can in queued_smiles or can == current_smiles
            needs_scheduling = can not in self._seen or (
                node is not None
                and not already_expanded
                and not already_scheduled
            )
            if needs_scheduling:
                pending.append(can)
                pending_set.add(can)
        return pending

    def is_complete(self) -> bool:
        return len(self._queue) == 0 and self._current_context is None

    def _unresolved_skips(self) -> List[Dict[str, str]]:
        unresolved: Dict[str, Dict[str, str]] = {}
        resolving_actions = {
            "accept-terminal",
            "commit",
        }
        for record in self.audit_state.decision_history:
            molecule = str(record.molecule or "").strip()
            action = str(record.action or "").strip()
            if not molecule:
                continue
            if action == "skip":
                unresolved[molecule] = {
                    "smiles": molecule,
                    "reason": str(
                        record.reasoning_detail or record.reasoning_summary or ""
                    ).strip(),
                }
            elif action in resolving_actions:
                unresolved.pop(molecule, None)
        return [unresolved[key] for key in sorted(unresolved)]

    def completion_blockers(self) -> Dict[str, Any]:
        """Return unresolved route state that forbids explicit finalization."""
        blockers: Dict[str, Any] = {}
        if self._current_context is not None:
            blockers["active_context"] = {
                "node_id": self._current_context.node_id,
                "smiles": self._current_context.smiles,
            }
        if self._queue:
            blockers["pending_queue_count"] = len(self._queue)
        if not self.tree.is_complete():
            blockers["incomplete_leaf_nodes"] = [
                {"node_id": node.node_id, "smiles": node.smiles}
                for node in self.tree.get_leaf_molecule_nodes()
                if node.role != MoleculeRole.TERMINAL.value
            ]
        unresolved_skips = self._unresolved_skips()
        if unresolved_skips:
            blockers["unresolved_skips"] = unresolved_skips
        return blockers

    def get_status(self) -> Dict[str, Any]:
        return {
            "target": self.tree.target,
            "status": self.tree.status,
            "steps_executed": self._steps_executed,
            "pending_count": self.pending_count,
            "total_molecules": len(self.tree.molecule_nodes),
            "total_steps": self.tree.total_steps,
            "max_depth": self.tree.max_depth,
            "elapsed_sec": round(time.time() - self._start_time, 1),
            "tree_complete": self.tree.is_complete(),
        }

    def peek_queue(self, n: int = 5) -> List[Dict[str, Any]]:
        """预览队列中前 n 个待处理分子。"""
        previews = []
        for i, (smiles, depth) in enumerate(self._queue):
            if i >= n:
                break
            node = self.tree.get_molecule_by_smiles(smiles)
            previews.append({
                "index": i,
                "smiles": smiles,
                "depth": depth,
                "cs_score": node.cs_score if node else 0,
            })
        return previews

    def select_next(self, smiles: str) -> bool:
        """将指定分子移到队列头部。"""
        can = canonical(smiles) or smiles
        for i, (smi, depth) in enumerate(self._queue):
            if smi == can:
                if i == 0:
                    return True
                del self._queue[i]
                self._queue.appendleft((can, depth))
                return True
        return False

    # ── 终止判定 ──

    def _check_terminal(self, smiles: str) -> Tuple[bool, Dict[str, Any]]:
        """三层终止判定。返回 (is_terminal, cs_result)。"""
        cs_result = compute_cs_score(
            smiles, knowledge_profile=self.knowledge_profile
        )

        # 层级 1: 极小分子
        mol = parse_mol(smiles)
        if mol and mol.GetNumHeavyAtoms() <= 6:
            cs_result["_terminal_reason"] = "small_molecule"
            return True, cs_result

        # 层级 2: CS score 阈值
        cs_score = cs_result.get("cs_score", 0)
        if (
            cs_score <= self.terminal_cs_threshold
            and not has_assigned_stereocenter(smiles)
        ):
            cs_result["_terminal_reason"] = "cs_threshold"
            return True, cs_result

        # 层级 3: 无可断键位（在 prepare_next 中检查）
        return False, cs_result


    # ── 主流程: prepare_next → LLM 决策 → commit_decision ──

    def prepare_next(self) -> Optional[ProposalContext]:
        """从队列取下一个分子，分析并返回 context 供 LLM 决策。"""
        expanded_products = {rxn.product_node for rxn in self.tree.reaction_nodes}

        # 如果上一个 context 未完成，放回队列
        if self._current_context is not None:
            old = self._current_context
            node = self.tree.get_molecule_by_smiles(old.smiles)
            if (
                node
                and node.role != MoleculeRole.TERMINAL.value
                and node.node_id not in expanded_products
            ):
                self._queue.appendleft((old.smiles, old.depth))
            self._current_context = None

        # max_steps 安全阀
        if self.step_budget_exhausted():
            logger.warning("max_steps=%d reached; unresolved queue retained", self.max_steps)
            return None

        while self._queue:
            smiles, depth = self._queue.popleft()
            node = self.tree.get_molecule_by_smiles(smiles)
            if node is None:
                break
            if node.role == MoleculeRole.TERMINAL.value:
                continue
            if node.node_id in expanded_products:
                continue
            break
        else:
            return None

        force_standard = smiles in self._force_standard_smiles

        # 终止判定
        is_terminal, cs_result = self._check_terminal(smiles)
        self.tree.update_complexity(smiles, cs_result)

        node = self.tree.get_molecule_by_smiles(smiles)
        is_target = (node.role == MoleculeRole.TARGET.value) if node else False

        # 快速通道: terminal
        if is_terminal and not is_target and not force_standard:
            self.tree.mark_terminal(smiles)
            ctx = ProposalContext(
                smiles=smiles,
                node_id=node.node_id if node else "",
                depth=depth,
                cs_score=cs_result.get("cs_score", 0),
                classification=cs_result.get("classification", ""),
                is_terminal=True,
                is_target=False,
                steps_executed=self._steps_executed,
                steps_remaining=max(0, self.max_steps - self._steps_executed),
                decision_tier="quick_pass",
                knowledge_profile=self.knowledge_profile,
            )
            self._current_context = ctx
            return ctx

        # 标准流程: 生成决策上下文
        decision_ctx = build_decision_context(
            smiles,
            knowledge_profile=self.knowledge_profile,
        )
        self.tree.update_decision_context(smiles, decision_ctx)

        # 层级 3: 无可断键位 → terminal
        bonds = decision_ctx.get("disconnectable_bonds", [])
        has_viable = any(b.get("alternatives") for b in bonds)
        fgi_options = decision_ctx.get("fgi_options", [])
        if not has_viable and not fgi_options:
            if not is_target and not force_standard:
                self.tree.mark_terminal(smiles)
                ctx = ProposalContext(
                    smiles=smiles,
                    node_id=node.node_id if node else "",
                    depth=depth,
                    cs_score=cs_result.get("cs_score", 0),
                    classification=cs_result.get("classification", ""),
                    is_terminal=True,
                    steps_executed=self._steps_executed,
                    steps_remaining=max(0, self.max_steps - self._steps_executed),
                    decision_tier="quick_pass",
                    knowledge_profile=self.knowledge_profile,
                )
                self._current_context = ctx
                return ctx

        # 深度限制标记（建议性）
        depth_limited = depth > self.max_depth

        ctx = ProposalContext(
            smiles=smiles,
            node_id=node.node_id if node else "",
            depth=depth,
            cs_score=cs_result.get("cs_score", 0),
            classification=cs_result.get("classification", ""),
            is_terminal=is_terminal,
            is_target=is_target,
            depth_limited=depth_limited,
            decision_context=decision_ctx,
            seen_smiles=set(self._seen),
            steps_executed=self._steps_executed,
            steps_remaining=max(0, self.max_steps - self._steps_executed),
            queue_preview=self.peek_queue(5),
            audit_state_summary=self.audit_state.get_summary(),
            failed_attempts_for_current=self.audit_state.get_failures_for_molecule(smiles),
            decision_tier="standard",
            knowledge_profile=self.knowledge_profile,
        )
        self._current_context = ctx
        return ctx

    def commit_decision(
        self,
        *,
        bond: Optional[Tuple[int, int]] = None,
        reaction_type: str = "",
        template_id: str = "",
        template_name: str = "",
        fgi_template_id: Optional[str] = None,
        fgi_template_name: str = "",
        precursor_override: Optional[List[str]] = None,
        reagents: Optional[List[str]] = None,
        validated_forward_validation: Optional[Dict[str, Any]] = None,
        llm_decision: Optional[LLMDecision] = None,
        llm_analysis: Optional[Dict[str, Any]] = None,
    ) -> CommitResult:
        """LLM 做完决策后，执行断键/FGI 并写入树。

        三种提交方式（优先级）:
          1. precursor_override — LLM 直接指定前体 SMILES
          2. fgi_template_id — FGI 操作
          3. bond + reaction_type — 断键操作
        """
        ctx = self._current_context
        if ctx is None:
            return CommitResult(success=False, error="no active context, call prepare_next first")

        smiles = ctx.smiles
        depth = ctx.depth

        # ── 执行断键/FGI ──
        precursors: List[str] = []
        actual_reaction_type = reaction_type
        actual_template_id = template_id

        if precursor_override is not None:
            # 方式 1: LLM 直接指定前体
            precursors = parse_precursors(precursor_override)
            if not precursors:
                return CommitResult(success=False, error="empty precursor_override")
            actual_reaction_type = reaction_type or "llm_proposed"

        elif fgi_template_id:
            # 方式 2: FGI
            fgi_result = execute_fgi(
                smiles,
                fgi_template_id,
                knowledge_profile=self.knowledge_profile,
            )
            if not fgi_result.success:
                self.audit_state.record_failure(
                    smiles, fgi_template_name or fgi_template_id,
                    fgi_result.error or "FGI failed",
                )
                self._current_context = None
                return CommitResult(success=False, error=fgi_result.error or "FGI failed")
            precursors = fgi_result.precursors
            actual_reaction_type = "fgi"
            actual_template_id = fgi_template_id

        elif bond is not None:
            # 方式 3: 断键
            break_result = None
            if template_id:
                mol = parse_mol(smiles)
                if mol is not None:
                    break_result = try_retro_template(
                        mol,
                        bond,
                        template_id,
                        knowledge_profile=self.knowledge_profile,
                    )
                if break_result is None or not break_result.success:
                    self.audit_state.record_failure(
                        smiles, template_name or template_id,
                        "exact template execution failed",
                    )
                    self._current_context = None
                    return CommitResult(success=False, error="exact template execution failed")
            else:
                break_result = execute_disconnection(
                    smiles,
                    bond,
                    reaction_type,
                    knowledge_profile=self.knowledge_profile,
                )
            if not break_result.success:
                self.audit_state.record_failure(
                    smiles, reaction_type,
                    break_result.error or "disconnection failed",
                )
                self._current_context = None
                return CommitResult(success=False, error=break_result.error or "disconnection failed")
            precursors = break_result.precursors

        else:
            self._current_context = None
            return CommitResult(success=False, error="must specify bond, fgi_template_id, or precursor_override")

        if not precursors:
            self._current_context = None
            return CommitResult(success=False, error="no precursors generated")

        preflight = _preflight_precursors(precursors)
        precursor_normalization = _precursor_normalization_payload(precursors, preflight)
        invalid_precursors = preflight.get("invalid", [])
        if invalid_precursors:
            self._current_context = None
            detail = _invalid_precursor_detail(invalid_precursors)
            return CommitResult(
                success=False,
                error=f"invalid precursor SMILES after preflight: {detail}",
            )
        precursors = list(preflight.get("precursors", []) or [])
        pending_precursors = self._pending_precursors_for_queue(precursors)
        if (
            pending_precursors
            and len(self._queue) + len(pending_precursors) > self.max_queue_size
        ):
            logger.warning(
                "Queue capacity exceeded: pending=%d new=%d max=%d",
                len(self._queue),
                len(pending_precursors),
                self.max_queue_size,
            )
            return CommitResult(success=False, error="queue capacity exceeded")

        # Record analysis only after every pre-mutation capacity check passes.
        if llm_analysis:
            self.tree.update_llm_analysis(smiles, llm_analysis)

        # ── 正向验证 ──
        fv_result = validated_forward_validation
        if fv_result is None and self.auto_forward_validate:
            try:
                fv_result = validate_forward(
                    precursors, smiles,
                    template_id=actual_template_id or None,
                    reaction_category=actual_reaction_type or None,
                    knowledge_profile=self.knowledge_profile,
                )
            except Exception as e:
                logger.warning("Forward validation error: %s", e)

        # ── 环路检测 ──
        cycle_warnings = []
        for smi in precursors:
            can = canonical(smi) or smi
            if can in self._seen:
                cycle_warnings.append(
                    f"Precursor {can} is already in the tree and may create a cycle."
                )

        # ── 写入树 ──
        evidence_refs = _collect_knowledge_refs(fv_result or {})
        if actual_template_id:
            try:
                template_ref = self.knowledge_profile.source(
                    "chem.reactions", actual_template_id
                )
            except Exception:
                template_ref = None
            if template_ref and template_ref not in evidence_refs:
                evidence_refs.insert(0, template_ref)

        evidence = TemplateEvidence(
            template_id=actual_template_id,
            template_name=template_name or fgi_template_name,
            reaction_category="fgi" if fgi_template_id is not None else actual_reaction_type,
            bond_atoms=list(bond) if bond else [],
            is_fgi=fgi_template_id is not None,
            knowledge_profile_hash=self.knowledge_profile.digest,
            knowledge_refs=evidence_refs,
        )

        rxn = self.tree.add_reaction(
            product_smiles=smiles,
            precursors=precursors,
            reagents=list(reagents or []),
            reaction_type=actual_reaction_type,
            template_evidence=evidence,
            llm_decision=llm_decision,
            forward_validation=fv_result,
            depth=depth,
        )

        self._steps_executed += 1
        self.audit_state.linear_step_count = self._steps_executed

        # 记录决策
        self.audit_state.record_decision(
            step_id=rxn.step_id,
            molecule=smiles,
            action="commit",
            reaction_name=template_name or fgi_template_name or actual_reaction_type,
            reasoning_summary=llm_decision.selection_reasoning if llm_decision else "",
            outcome="committed",
            confidence=llm_decision.confidence if llm_decision else "medium",
        )

        # ── 分类前体 ──
        new_pending: List[str] = []
        new_terminal: List[str] = []
        pending_precursor_set = set(pending_precursors)

        for smi in precursors:
            can = canonical(smi) or smi
            is_term, cs_result = self._check_terminal(can)
            self.tree.update_complexity(can, cs_result)

            if is_term:
                self.tree.mark_terminal(can)
                new_terminal.append(can)
            else:
                self.tree.mark_intermediate(can)
                node = self.tree.get_molecule_by_smiles(can)
                if can in pending_precursor_set:
                    self._seen.add(can)
                    queue_depth = node.depth if node is not None else depth + 1
                    self._queue.append((can, queue_depth))
                    new_pending.append(can)
                else:
                    logger.info("Convergent or already scheduled: %s", can)

        self._current_context = None

        return CommitResult(
            success=True,
            reaction_node=rxn,
            new_pending=new_pending,
            new_terminal=new_terminal,
            tree_complete=self.is_complete() and self.tree.is_complete(),
            cycle_warnings=cycle_warnings,
            forward_validation=fv_result,
            knowledge_profile_hash=self.knowledge_profile.digest,
            knowledge_refs=evidence_refs,
        )


    # ── 便捷方法 ──

    def accept_terminal(self, smiles: Optional[str] = None, reason: str = "") -> None:
        """LLM 确认标记分子为 terminal。"""
        target = smiles
        if target is None and self._current_context:
            target = self._current_context.smiles

        if not target:
            return

        can = canonical(target) or target
        self.tree.mark_terminal(can)

        self.audit_state.record_decision(
            step_id="",
            molecule=can,
            action="accept-terminal",
            reasoning_summary=reason[:120],
            reasoning_detail=reason,
            outcome="terminal",
        )

        # 清除 context（如果是当前分子）
        if self._current_context and self._current_context.smiles == can:
            self._current_context = None

    def skip_current(self, reason: str = "no viable proposals") -> None:
        """跳过当前分子，标记为 terminal。"""
        if not self._current_context:
            return

        smiles = self._current_context.smiles
        self.tree.mark_terminal(smiles)

        self.audit_state.record_decision(
            step_id="",
            molecule=smiles,
            action="skip",
            reasoning_summary=reason[:120],
            reasoning_detail=reason,
            outcome="skipped",
        )
        self._current_context = None

    # ── 分层上下文: 按需展开 ──

    def explore_bond(self, bond_idx: int) -> Dict[str, Any]:
        """按需展开某个键位的完整前体方案。

        LLM 在 compact 视图中看到键位概览后，对感兴趣的键位调用此方法
        获取完整的前体 SMILES 和模板信息。
        """
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context"}

        bonds = ctx.decision_context.get("disconnectable_bonds", [])
        if bond_idx < 0 or bond_idx >= len(bonds):
            return {"error": f"bond_idx {bond_idx} out of range (0-{len(bonds)-1})"}

        b = bonds[bond_idx]
        # Keep the original detail payload as a comment instead of deleting it.
        # Reason: we are adding only the two lightweight audit identifiers that
        # help verify site identity, while leaving the rest of explore_bond intact.
        # result = {
        #     "bond_idx": bond_idx,
        #     "atoms": b.get("atoms", []),
        #     "bond_type": b.get("bond_type", ""),
        #     "in_ring": b.get("in_ring", False),
        #     "heuristic_score": b.get("heuristic_score", 0),
        #     "alternatives": b.get("alternatives", []),
        #     "hint": "用 try_disconnection(bond_idx, alt_idx) 沙盒试断。",
        # }
        bond_fg_context = _build_bond_fg_context(ctx.smiles, tuple(b.get("atoms", [])))
        result = {
            "bond_idx": bond_idx,
            "actual_bond_idx": b.get("actual_bond_idx", -1),
            "atoms": b.get("atoms", []),
            "role_pair": b.get("role_pair", ["unclear", "unclear"]),
            # New minimal local audit: keep full molecular maps out of compact
            # decision_context, but expose the bond-specific FG micro-context here
            # where the model is actively checking site identity and mechanism fit.
            "bond_fg_context": bond_fg_context,
            "bond_type": b.get("bond_type", ""),
            "in_ring": b.get("in_ring", False),
            "heuristic_score": b.get("heuristic_score", 0),
            "alternatives": b.get("alternatives", []),
            "hint": "Use try_action to test this disconnection in the sandbox.",
        }
        # 包含 smart_capping（如果 build_decision_context 已生成）
        if "smart_capping" in b:
            result["smart_capping"] = b["smart_capping"]
        else:
            # 按需生成
            try:
                from Rachel.chem_tools.smart_cap import suggest_capping
                cap_result = suggest_capping(
                    ctx.smiles,
                    tuple(b.get("atoms", [])),
                    knowledge_profile=self.knowledge_profile,
                )
                if cap_result.get("ok") and cap_result.get("proposals"):
                    result["smart_capping"] = cap_result["proposals"][:3]
            except Exception:
                pass
        return result

    def explore_fgi(self) -> Dict[str, Any]:
        """按需展开所有 FGI 选项的完整信息。"""
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context"}

        fgi_list = ctx.decision_context.get("fgi_options", [])
        return {
            "n_fgi": len(fgi_list),
            "fgi_options": fgi_list,
            "hint": "Use try_action to test this FGI in the sandbox.",
        }

    # ── 沙盒试探: 不写入树 ──

    def explore_reaction(self, reaction_id: str) -> Dict[str, Any]:
        """Expand one concrete reaction family while preserving candidate IDs."""
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context"}
        result = expand_reaction_family(ctx.decision_context, reaction_id)
        result["prompt_brief"] = build_prompt_brief(
            build_prompt_mount(
                "explore_reaction",
                payload=result,
                knowledge_profile=self.knowledge_profile,
            )
        )
        return result

    def reaction_sites(self) -> Dict[str, Any]:
        """Return the complete first-layer site/reaction menu."""
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context"}
        site_map = build_site_reaction_map(ctx.decision_context)
        return {
            "primary_disclosure": "site_reaction_map",
            "site_expand_command": "explore_site(site_id)",
            "action_count_meaning": "Number of concrete action instances at this site.",
            "site_reaction_map": site_map,
            "site_count": len(site_map),
            "total_reaction_count": sum(
                len(site.get("reactions", []) or []) for site in site_map
            ),
            "next_step": "explore_site(site_id)",
            "hint": "Choose a site_id, then call explore_site(site_id).",
            "prompt_brief": build_prompt_brief(
                build_prompt_mount(
                    "reaction_sites",
                    decision_context=ctx.decision_context,
                    payload={"site_reaction_map": site_map},
                    knowledge_profile=self.knowledge_profile,
                )
            ),
        }

    def explore_site(self, site_id: str) -> Dict[str, Any]:
        """Expand all normalized candidates that compete at one real reaction site."""
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context"}
        result = expand_site_candidates(ctx.decision_context, site_id)
        if result.get("site_type") == "bond":
            atoms = result.get("atoms", [])
            if isinstance(atoms, list) and len(atoms) == 2:
                result["bond_fg_context"] = _build_bond_fg_context(ctx.smiles, tuple(atoms))
        result["prompt_brief"] = build_prompt_brief(
            build_prompt_mount(
                "explore_site",
                payload=result,
                knowledge_profile=self.knowledge_profile,
            )
        )
        return result

    def _sandbox_classify_precursors(self, precursors: List[str]) -> Tuple[
        List[Dict[str, Any]], List[str]
    ]:
        """对前体做 CS 评分和 terminal 判定（不写入树）。"""
        details = []
        cycle_warnings = []
        for smi in precursors:
            can = canonical(smi) or smi
            cs = compute_cs_score(
                can, knowledge_profile=self.knowledge_profile
            )
            is_term, _ = self._check_terminal(can)
            mol = parse_mol(can)
            entry: Dict[str, Any] = {
                "smiles": can,
                "cs_score": cs.get("cs_score", 0),
                "classification": cs.get("classification", ""),
                "is_terminal": is_term,
                "heavy_atoms": mol.GetNumHeavyAtoms() if mol else 0,
            }
            if can in self._seen:
                cycle_warnings.append(f"Precursor {can} is already in the tree.")
                entry["already_seen"] = True
            details.append(entry)
        return details, cycle_warnings

    def try_disconnection(
        self,
        bond_idx: int = 0,
        alt_idx: int = 0,
        *,
        bond: Optional[Tuple[int, int]] = None,
        reaction_type: str = "",
        template_id: str = "",
        allow_template_fallback: bool = True,
        action_context: Optional[Dict[str, Any]] = None,
    ) -> SandboxResult:
        """沙盒试断键。不写入树，返回前体 + 验证结果。

        两种调用方式:
          1. bond_idx + alt_idx — 从 explore_bond 的结果中选
          2. bond + reaction_type — 直接指定原子对和反应类型
        """
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return SandboxResult(success=False, error="no active context")

        smiles = ctx.smiles

        # 解析参数
        actual_bond = bond
        actual_type = reaction_type
        actual_tid = ""
        actual_tname = ""
        execution_mode = ""
        executed_template_id = ""
        bond_meta: Dict[str, Any] = {}
        alt_meta: Dict[str, Any] = {}

        if bond is None:
            bonds = ctx.decision_context.get("disconnectable_bonds", [])
            if bond_idx < 0 or bond_idx >= len(bonds):
                return SandboxResult(
                    success=False,
                    error=f"bond_idx {bond_idx} out of range",
            )
            b = bonds[bond_idx]
            bond_meta = dict(b)
            actual_bond = tuple(b["atoms"])
            alts = b.get("alternatives", [])
            if alt_idx < 0 or alt_idx >= len(alts):
                return SandboxResult(
                    success=False,
                    error=f"alt_idx {alt_idx} out of range (0-{len(alts)-1})",
            )
            alt = alts[alt_idx]
            alt_meta = dict(alt)
            actual_type = alt.get("template", "").split("(")[0].strip()
            actual_tid = alt.get("template_id", "")
            actual_tname = alt.get("template", "")
        if template_id:
            actual_tid = template_id

        # 执行断键
        break_result = None
        if actual_tid:
            mol = parse_mol(smiles)
            if mol is not None:
                break_result = try_retro_template(
                    mol,
                    actual_bond,
                    actual_tid,
                    knowledge_profile=self.knowledge_profile,
                )
            if break_result is not None and break_result.success:
                execution_mode = "exact_template"
                executed_template_id = actual_tid
            elif not allow_template_fallback:
                return SandboxResult(
                    success=False,
                    error=f"exact template execution failed: {actual_tid}",
                    reaction_type=actual_type,
                    template_id=actual_tid,
                    template_name=actual_tname,
                    execution_mode="exact_template_failed",
                    declared_template_id=actual_tid,
                    candidate_consistency={
                        "declared_template_id": actual_tid,
                        "preview_matches_sandbox": False,
                        "fallback_used": False,
                    },
                )

        if break_result is None or not break_result.success:
            break_result = execute_disconnection(
                smiles,
                actual_bond,
                actual_type,
                knowledge_profile=self.knowledge_profile,
            )
            execution_mode = "fuzzy_fallback" if actual_tid else "reaction_type"
        if not break_result.success:
            return SandboxResult(
                success=False,
                error=break_result.error or "disconnection failed",
                reaction_type=actual_type,
                template_id=actual_tid,
                template_name=actual_tname,
                execution_mode=execution_mode,
                declared_template_id=actual_tid,
            )

        precursors = break_result.precursors
        if not precursors:
            return SandboxResult(success=False, error="no precursors generated")

        original_precursors = list(precursors)
        preflight = _preflight_precursors(original_precursors)
        precursor_normalization = _precursor_normalization_payload(original_precursors, preflight)
        invalid_precursors = preflight.get("invalid", [])
        if invalid_precursors:
            detail = _invalid_precursor_detail(invalid_precursors)
            return SandboxResult(
                success=False,
                error=f"invalid precursor SMILES after preflight: {detail}",
                precursors=original_precursors,
                reaction_type=actual_type,
                template_id=actual_tid,
                template_name=actual_tname,
                execution_mode=execution_mode,
                declared_template_id=actual_tid,
                executed_template_id=executed_template_id,
                candidate_consistency={
                    "declared_template_id": actual_tid,
                    "executed_template_id": executed_template_id,
                    "preview_matches_sandbox": execution_mode == "exact_template",
                    "fallback_used": execution_mode == "fuzzy_fallback",
                } if actual_tid else None,
                precursor_normalization=precursor_normalization,
            )
        precursors = list(preflight.get("precursors", []) or [])

        # 前体分析
        details, cycle_warnings = self._sandbox_classify_precursors(precursors)

        # 正向验证
        fv = None
        validation_action_context = dict(action_context or {})
        validation_action_context.setdefault("source", "bond")
        validation_action_context.setdefault("template_id", actual_tid)
        validation_action_context.setdefault("template_name", actual_tname)
        validation_action_context.setdefault("reaction_type", actual_type)
        validation_action_context.setdefault("actual_bond_idx", bond_meta.get("actual_bond_idx"))
        validation_action_context.setdefault("atoms", list(actual_bond or []))
        validation_action_context.setdefault("bond_type", bond_meta.get("bond_type", ""))
        validation_action_context.setdefault("in_ring", bool(bond_meta.get("in_ring", False)))
        risk_tags = list(validation_action_context.get("risk_tags", []) or [])
        if validation_action_context.get("in_ring") and "ring_bond" not in risk_tags:
            risk_tags.append("ring_bond")
        if alt_meta.get("incompatible_with") and "functional_group_compatibility" not in risk_tags:
            risk_tags.append("functional_group_compatibility")
        if precursor_normalization and "organometallic_source_obligation" not in risk_tags:
            risk_tags.append("organometallic_source_obligation")
        validation_action_context["risk_tags"] = risk_tags
        if precursor_normalization:
            validation_action_context["precursor_normalization"] = precursor_normalization
        try:
            fv = validate_forward(
                precursors, smiles,
                template_id=actual_tid or None,
                reaction_category=actual_type or None,
                action_context=validation_action_context,
                knowledge_profile=self.knowledge_profile,
            )
        except Exception as exc:
            fv = _validation_unavailable_payload(
                exc,
                precursors,
                smiles,
                validation_action_context,
            )

        # 原子平衡
        ab = None
        try:
            ab = check_atom_balance(precursors, smiles)
        except Exception:
            pass

        return SandboxResult(
            success=True,
            precursors=precursors,
            precursor_details=details,
            forward_validation=fv,
            validation_micro=_build_validation_micro(fv),
            atom_balance=ab,
            cycle_warnings=cycle_warnings,
            reaction_type=actual_type,
            eas_site_audit=(fv.get("checks") or {}).get("eas_site_selectivity"),
            template_id=actual_tid,
            template_name=actual_tname,
            execution_mode=execution_mode,
            declared_template_id=actual_tid,
            executed_template_id=executed_template_id,
            precursor_normalization=precursor_normalization,
            candidate_consistency={
                "declared_template_id": actual_tid,
                "executed_template_id": executed_template_id,
                "preview_matches_sandbox": execution_mode == "exact_template",
                "fallback_used": execution_mode == "fuzzy_fallback",
            } if actual_tid else None,
        )

    def try_fgi(
        self,
        fgi_idx: int = 0,
        action_context: Optional[Dict[str, Any]] = None,
    ) -> SandboxResult:
        """沙盒试 FGI。不写入树。"""
        ctx = self._current_context
        if ctx is None or ctx.decision_context is None:
            return SandboxResult(success=False, error="no active context")

        smiles = ctx.smiles
        fgi_list = ctx.decision_context.get("fgi_options", [])
        if fgi_idx < 0 or fgi_idx >= len(fgi_list):
            return SandboxResult(
                success=False,
                error=f"fgi_idx {fgi_idx} out of range (0-{len(fgi_list)-1})",
            )

        fgi = fgi_list[fgi_idx]
        tid = fgi["template_id"]
        tname = fgi.get("template", "")

        fgi_result = execute_fgi(
            smiles,
            tid,
            knowledge_profile=self.knowledge_profile,
        )
        if not fgi_result.success:
            return SandboxResult(
                success=False,
                error=fgi_result.error or "FGI failed",
                template_name=tname,
            )

        precursors = fgi_result.precursors
        preflight = _preflight_precursors(precursors)
        precursor_normalization = _precursor_normalization_payload(precursors, preflight)
        invalid_precursors = preflight.get("invalid", [])
        if invalid_precursors:
            detail = _invalid_precursor_detail(invalid_precursors)
            return SandboxResult(
                success=False,
                error=f"invalid precursor SMILES after preflight: {detail}",
                precursors=precursors,
                reaction_type="fgi",
                template_id=tid,
                template_name=tname,
                precursor_normalization=precursor_normalization,
            )
        precursors = list(preflight.get("precursors", []) or [])

        details, cycle_warnings = self._sandbox_classify_precursors(precursors)

        fv = None
        validation_action_context = dict(action_context or {})
        validation_action_context.setdefault("source", "fgi")
        validation_action_context.setdefault("template_id", tid)
        validation_action_context.setdefault("template_name", tname)
        validation_action_context.setdefault("reaction_type", "fgi")
        validation_action_context.setdefault("in_ring", False)
        try:
            fv = validate_forward(
                precursors,
                smiles,
                template_id=tid,
                action_context=validation_action_context,
                knowledge_profile=self.knowledge_profile,
            )
        except Exception as exc:
            fv = _validation_unavailable_payload(
                exc,
                precursors,
                smiles,
                validation_action_context,
            )

        return SandboxResult(
            success=True,
            precursors=precursors,
            precursor_details=details,
            forward_validation=fv,
            validation_micro=_build_validation_micro(fv),
            cycle_warnings=cycle_warnings,
            reaction_type="fgi",
            eas_site_audit=(fv.get("checks") or {}).get("eas_site_selectivity"),
            template_id=tid,
            template_name=tname,
            precursor_normalization=precursor_normalization,
        )

    def try_precursors(
        self,
        precursors: List[str],
        reaction_type: str = "llm_proposed",
        action_context: Optional[Dict[str, Any]] = None,
        reagents: Optional[List[str]] = None,
    ) -> SandboxResult:
        """LLM 自己提出前体 SMILES，沙盒验证。不写入树。

        LLM 可以基于化学知识直接提出前体，编排器帮它验证:
          - SMILES 合法性
          - CS 评分 + terminal 判定
          - 正向验证（前体能否合成产物）
          - 原子平衡
          - 环路检测
        """
        ctx = self._current_context
        if ctx is None:
            return SandboxResult(success=False, error="no active context")

        smiles = ctx.smiles

        # 解析前体
        parsed = parse_precursors(precursors)
        if not parsed:
            return SandboxResult(success=False, error="empty or invalid precursors")

        # 验证 SMILES 合法性
        preflight = _preflight_precursors(parsed)
        precursor_normalization = _precursor_normalization_payload(parsed, preflight)
        invalid_precursors = preflight.get("invalid", [])
        if invalid_precursors:
            detail = _invalid_precursor_detail(invalid_precursors)
            return SandboxResult(
                success=False,
                error=f"invalid precursor SMILES after preflight: {detail}",
                precursor_normalization=precursor_normalization,
            )
        valid_precursors = list(preflight.get("precursors", []) or [])

        reagent_preflight = _preflight_precursors(list(reagents or []))
        invalid_reagents = reagent_preflight.get("invalid", [])
        if invalid_reagents:
            detail = _invalid_precursor_detail(invalid_reagents)
            return SandboxResult(
                success=False,
                error=f"invalid reagent SMILES after preflight: {detail}",
            )
        valid_reagents = list(reagent_preflight.get("precursors", []) or [])

        # 前体分析
        details, cycle_warnings = self._sandbox_classify_precursors(valid_precursors)

        # 正向验证
        fv = None
        validation_action_context = dict(action_context or {})
        validation_action_context.setdefault("source", "custom_precursors")
        validation_action_context.setdefault("reaction_type", reaction_type)
        validation_action_context["reagents"] = valid_reagents
        if "in_ring" not in validation_action_context:
            validation_action_context["in_ring"] = _custom_context_implies_ring_action(
                validation_action_context,
                reaction_type,
            )
        if precursor_normalization:
            validation_action_context["precursor_normalization"] = precursor_normalization
            risk_tags = list(validation_action_context.get("risk_tags", []) or [])
            if "organometallic_source_obligation" not in risk_tags:
                risk_tags.append("organometallic_source_obligation")
            validation_action_context["risk_tags"] = risk_tags
        try:
            fv = validate_forward(
                valid_precursors, smiles,
                reaction_category=reaction_type,
                action_context=validation_action_context,
                reagents=valid_reagents,
                knowledge_profile=self.knowledge_profile,
            )
        except Exception as exc:
            fv = _validation_unavailable_payload(
                exc,
                valid_precursors,
                smiles,
                validation_action_context,
            )

        # 原子平衡
        ab = None
        try:
            ab = check_atom_balance(valid_precursors + valid_reagents, smiles)
        except Exception:
            pass

        site_audit = audit_site_retention(smiles, valid_precursors)

        return SandboxResult(
            success=True,
            precursors=valid_precursors,
            reagents=valid_reagents,
            precursor_details=details,
            forward_validation=fv,
            validation_micro=_build_validation_micro(fv),
            atom_balance=ab,
            cycle_warnings=cycle_warnings,
            reaction_type=reaction_type,
            site_audit=site_audit,
            eas_site_audit=(fv.get("checks") or {}).get("eas_site_selectivity"),
            precursor_normalization=precursor_normalization,
        )

    def finalize(self, llm_summary: str = "") -> Dict[str, Any]:
        """完成编排，返回最终 JSON。"""
        blockers = self.completion_blockers()
        if blockers:
            return {
                "error": "route_not_ready_for_finalize",
                "completion_blockers": blockers,
            }

        self.tree.complete(llm_summary)

        return {
            "status": self.get_status(),
            "tree": self.tree.to_dict(),
            "audit_state": self.audit_state.to_dict(),
        }

    def get_tree(self) -> RetrosynthesisTree:
        return self.tree

    # ── 自动模式（测试用） ──

    def auto_run(self, verbose: bool = True) -> Dict[str, Any]:
        """全自动运行（贪心策略，不需要 LLM）。用于测试。"""
        iteration = 0
        max_iter = self.max_steps * 2

        while iteration < max_iter:
            iteration += 1

            ctx = self.prepare_next()
            if ctx is None:
                break

            # 快速通道
            if ctx.decision_tier == "quick_pass":
                if verbose:
                    print(f"  ✓ terminal: {ctx.smiles[:50]}  CS={ctx.cs_score:.2f}")
                self.accept_terminal(reason="auto: quick_pass")
                continue

            # 标准流程: 自动选最佳方案
            dc = ctx.decision_context
            if not dc:
                self.skip_current("no decision context")
                continue

            if verbose:
                print(f"\n[Step {self._steps_executed + 1}] depth={ctx.depth}  "
                      f"{ctx.smiles[:60]}  CS={ctx.cs_score:.2f}")

            # 尝试断键
            committed = False
            bonds = dc.get("disconnectable_bonds", [])
            sorted_bonds = sorted(
                bonds, key=lambda b: b.get("heuristic_score", 0), reverse=True,
            )
            for bond_info in sorted_bonds:
                alts = bond_info.get("alternatives", [])
                if not alts:
                    continue
                best_alt = alts[0]
                result = self.commit_decision(
                    bond=tuple(bond_info["atoms"]),
                    reaction_type=best_alt.get("template", "").split("(")[0].strip(),
                    template_id=best_alt.get("template_id", ""),
                    template_name=best_alt.get("template", ""),
                )
                if result.success:
                    committed = True
                    if verbose:
                        prec = " + ".join(
                            result.reaction_node.reaction_smiles.split(">>")[0].split(".")
                        ) if result.reaction_node else ""
                        print(f"  → {result.reaction_node.reaction_type}: {prec[:70]}")
                        for t in result.new_terminal:
                            print(f"    ✓ {t[:50]}")
                        for p in result.new_pending:
                            print(f"    … {p[:50]}")
                    break

            if not committed:
                # 尝试 FGI
                for fgi in dc.get("fgi_options", []):
                    result = self.commit_decision(
                        fgi_template_id=fgi["template_id"],
                        fgi_template_name=fgi.get("template", ""),
                    )
                    if result.success:
                        committed = True
                        if verbose:
                            print(f"  → FGI: {fgi.get('template', '')}")
                        break

            if not committed:
                if verbose:
                    print(f"  ✗ no viable disconnection")
                self.skip_current("auto: no viable disconnection")

        # 完成
        report = self.finalize("auto_run completed")

        if verbose:
            print(f"\n{'='*60}")
            status = self.get_status()
            print(f"Complete: steps={status['steps_executed']}  "
                  f"molecules={status['total_molecules']}  "
                  f"depth={status['max_depth']}  "
                  f"elapsed={status['elapsed_sec']}s")

        return report
