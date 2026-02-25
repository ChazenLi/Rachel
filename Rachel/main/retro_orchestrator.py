"""
BFS 逆合成编排引擎
==================
流程控制器，不做化学决策。所有判断信息传给 LLM，由 LLM 决定。

核心循环（分层上下文 + 沙盒试探）:
    prepare_next()       → 精简摘要（top-N 方案概览，不含完整前体）
    explore_bond(idx)    → 按需展开某个键位的详细方案
    explore_fgi()        → 按需展开 FGI 选项
    try_disconnection()  → 沙盒试断（不写入树，返回前体 + 验证）
    try_fgi()            → 沙盒试 FGI
    try_precursors()     → LLM 自提前体，沙盒验证
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

from Rachel.chem_tools._rdkit_utils import canonical, parse_mol
from Rachel.chem_tools.cs_score import compute_cs_score, classify_complexity
from Rachel.chem_tools.bond_break import (
    execute_disconnection, execute_fgi, preview_disconnections,
)
from Rachel.chem_tools.forward_validate import validate_forward, check_atom_balance
from Rachel.chem_tools.fg_warnings import suggest_protection_needs

from Rachel.tools.llm_retro_platform import build_decision_context, format_context_text

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

logger = logging.getLogger(__name__)


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

    def to_dict(self, detail: str = "compact", top_n: int = 5) -> Dict[str, Any]:
        """分层输出。

        detail:
          "compact" — 默认，键位只给概览（atoms, score, reaction_types, n_alternatives）
          "full"    — 完整版，包含所有前体方案（兼容旧行为，慎用）
        top_n: compact 模式下只展示前 N 个键位
        """
        if self.decision_tier == "quick_pass":
            return {
                "action": "awaiting_decision",
                "decision_tier": "quick_pass",
                "smiles": self.smiles,
                "node_id": self.node_id,
                "depth": self.depth,
                "cs_score": self.cs_score,
                "classification": self.classification,
                "is_terminal": self.is_terminal,
                "steps_executed": self.steps_executed,
                "hint": "快速通道：直接 accept_terminal 即可。",
            }

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
            d["molecule"] = ctx.get("molecule", {})
            d["functional_groups"] = ctx.get("functional_groups", [])
            d["complexity"] = ctx.get("complexity", {})
            d["n_bonds"] = ctx.get("n_bonds", 0)
            d["n_fgi"] = len(ctx.get("fgi_options", []))
            d["warnings"] = ctx.get("warnings", [])

            bonds = ctx.get("disconnectable_bonds", [])
            if detail == "full":
                # 完整版：所有键位 + 所有前体方案
                d["disconnectable_bonds"] = bonds
                d["fgi_options"] = ctx.get("fgi_options", [])
            else:
                # compact：只给键位概览，不含前体 SMILES
                bond_summary = []
                for i, b in enumerate(bonds[:top_n]):
                    alts = b.get("alternatives", [])
                    summary = {
                        "bond_idx": i,
                        "atoms": b.get("atoms", []),
                        "bond_type": b.get("bond_type", ""),
                        "in_ring": b.get("in_ring", False),
                        "heuristic_score": b.get("heuristic_score", 0),
                        "n_alternatives": len(alts),
                        "reaction_types": [
                            a.get("template", "").split("(")[0].strip()
                            for a in alts[:3]
                        ],
                    }
                    # 包含 smart_capping 概览（如果有）
                    sc = b.get("smart_capping", [])
                    if sc:
                        summary["smart_capping"] = [
                            {"type": c["reaction_type"], "conf": c["confidence"]}
                            for c in sc[:2]
                        ]
                    bond_summary.append(summary)
                d["bond_summary"] = bond_summary
                if len(bonds) > top_n:
                    d["bonds_omitted"] = len(bonds) - top_n

                # FGI 也只给概览
                fgi_list = ctx.get("fgi_options", [])
                if fgi_list:
                    d["fgi_summary"] = [
                        {"fgi_idx": i, "template": f.get("template", "")}
                        for i, f in enumerate(fgi_list[:5])
                    ]
                    if len(fgi_list) > 5:
                        d["fgi_omitted"] = len(fgi_list) - 5

                d["hint"] = (
                    "这是精简视图。用 explore_bond(idx) 查看某键位的完整前体方案，"
                    "用 try_disconnection() 沙盒试断，满意后 commit_decision()。"
                )

        if self.queue_preview:
            d["queue_preview"] = self.queue_preview
        if self.audit_state_summary:
            d["audit_state_summary"] = self.audit_state_summary
        if self.failed_attempts_for_current:
            d["failed_attempts_for_current"] = self.failed_attempts_for_current
        return d

    def to_text(self, target_smiles: str = "", detail: str = "compact") -> str:
        """格式化为 LLM 可读文本。"""
        if self.decision_tier == "quick_pass":
            return (
                f"分子 {self.smiles} — CS={self.cs_score:.2f} "
                f"({self.classification})，可作为终端原料。"
            )
        parts = []
        parts.append(f"═══ 逆合成决策 — 第 {self.steps_executed + 1} 步 ═══")
        if target_smiles:
            parts.append(f"目标产物: {target_smiles}")
        parts.append(f"当前分子: {self.smiles}  depth={self.depth}")
        parts.append(f"复杂度: CS={self.cs_score:.2f} ({self.classification})")
        parts.append("")

        if detail == "full" and self.decision_context:
            parts.append(format_context_text(self.decision_context))
        elif self.decision_context:
            # compact 文本：只列键位概览
            ctx = self.decision_context
            fgs = ctx.get("functional_groups", [])
            if fgs:
                fg_names = [g["name"] for g in fgs[:10]]
                parts.append(f"官能团: {', '.join(fg_names)}")

            bonds = ctx.get("disconnectable_bonds", [])
            parts.append(f"\n可断键位: {len(bonds)} 个")
            for i, b in enumerate(bonds[:5]):
                alts = b.get("alternatives", [])
                types = [a.get("template", "").split("(")[0].strip() for a in alts[:3]]
                parts.append(
                    f"  [{i}] atoms={b['atoms']}  score={b.get('heuristic_score', 0):.3f}"
                    f"  方案={len(alts)}  ({', '.join(types)})"
                )
            if len(bonds) > 5:
                parts.append(f"  ... 还有 {len(bonds) - 5} 个键位")

            fgi = ctx.get("fgi_options", [])
            if fgi:
                parts.append(f"\nFGI 选项: {len(fgi)} 个")
                for i, f in enumerate(fgi[:3]):
                    parts.append(f"  [{i}] {f.get('template', '')}")

            parts.append("\n提示: explore_bond(idx) 查看详情, try_disconnection() 沙盒试断")

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
    precursor_details: List[Dict[str, Any]] = field(default_factory=list)
    forward_validation: Optional[Dict[str, Any]] = None
    atom_balance: Optional[Dict[str, Any]] = None
    cycle_warnings: List[str] = field(default_factory=list)
    reaction_type: str = ""
    template_id: str = ""
    template_name: str = ""
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {"success": self.success}
        if self.precursors:
            d["precursors"] = self.precursors
        if self.precursor_details:
            d["precursor_details"] = self.precursor_details
        if self.forward_validation:
            from .retro_tree import _flatten_fv
            d["forward_validation"] = _flatten_fv(self.forward_validation)
        if self.atom_balance:
            d["atom_balance"] = {
                "balanced": self.atom_balance.get("balanced", False),
                "message": self.atom_balance.get("message", ""),
            }
        if self.cycle_warnings:
            d["cycle_warnings"] = self.cycle_warnings
        if self.reaction_type:
            d["reaction_type"] = self.reaction_type
        if self.template_name:
            d["template_name"] = self.template_name
        if self.error:
            d["error"] = self.error
        d["hint"] = "这是沙盒结果，未写入树。满意请调 commit_decision()。"
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
        terminal_cs_threshold: float = 2.5,
        auto_forward_validate: bool = True,
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

        # BFS 队列: (smiles, depth)
        self._queue: deque = deque()
        self._queue.append((can, 0))

        # 已见分子集合
        self._seen: Set[str] = {can}

        # 当前活跃 context
        self._current_context: Optional[ProposalContext] = None

        # 统计
        self._steps_executed: int = 0
        self._start_time: float = time.time()

        # 审计状态
        self.audit_state = SynthesisAuditState()

        # 分析目标分子复杂度
        cs = compute_cs_score(can)
        self.tree.update_complexity(can, cs)
        self.audit_state.set_target_complexity(cs.get("cs_score", 0))

    # ── 状态查询 ──

    @property
    def pending_count(self) -> int:
        return len(self._queue)

    def is_complete(self) -> bool:
        return len(self._queue) == 0 and self._current_context is None

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
        cs_result = compute_cs_score(smiles)

        # 层级 1: 极小分子
        mol = parse_mol(smiles)
        if mol and mol.GetNumHeavyAtoms() <= 6:
            cs_result["_terminal_reason"] = "small_molecule"
            return True, cs_result

        # 层级 2: CS score 阈值
        cs_score = cs_result.get("cs_score", 0)
        if cs_score <= self.terminal_cs_threshold:
            cs_result["_terminal_reason"] = "cs_threshold"
            return True, cs_result

        # 层级 3: 无可断键位（在 prepare_next 中检查）
        return False, cs_result


    # ── 主流程: prepare_next → LLM 决策 → commit_decision ──

    def prepare_next(self) -> Optional[ProposalContext]:
        """从队列取下一个分子，分析并返回 context 供 LLM 决策。"""
        # 如果上一个 context 未完成，放回队列
        if self._current_context is not None:
            old = self._current_context
            node = self.tree.get_molecule_by_smiles(old.smiles)
            if node and node.role != MoleculeRole.TERMINAL.value:
                self._queue.appendleft((old.smiles, old.depth))
            self._current_context = None

        if not self._queue:
            return None

        # max_steps 安全阀
        if self._steps_executed >= self.max_steps:
            logger.warning("max_steps=%d reached, draining queue", self.max_steps)
            while self._queue:
                smi, _ = self._queue.popleft()
                self.tree.mark_terminal(smi)
            return None

        smiles, depth = self._queue.popleft()

        # 终止判定
        is_terminal, cs_result = self._check_terminal(smiles)
        self.tree.update_complexity(smiles, cs_result)

        node = self.tree.get_molecule_by_smiles(smiles)
        is_target = (node.role == MoleculeRole.TARGET.value) if node else False

        # 快速通道: terminal
        if is_terminal and not is_target:
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
            )
            self._current_context = ctx
            return ctx

        # 标准流程: 生成决策上下文
        decision_ctx = build_decision_context(smiles)
        self.tree.update_decision_context(smiles, decision_ctx)

        # 层级 3: 无可断键位 → terminal
        bonds = decision_ctx.get("disconnectable_bonds", [])
        has_viable = any(b.get("alternatives") for b in bonds)
        fgi_options = decision_ctx.get("fgi_options", [])
        if not has_viable and not fgi_options:
            if not is_target:
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

        # 记录 LLM 分析
        if llm_analysis:
            self.tree.update_llm_analysis(smiles, llm_analysis)

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
            fgi_result = execute_fgi(smiles, fgi_template_id)
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
            break_result = execute_disconnection(smiles, bond, reaction_type)
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

        # ── 正向验证 ──
        fv_result = None
        if self.auto_forward_validate:
            try:
                fv_result = validate_forward(
                    precursors, smiles,
                    template_id=actual_template_id or None,
                    reaction_category=actual_reaction_type or None,
                )
            except Exception as e:
                logger.warning("Forward validation error: %s", e)

        # ── 环路检测 ──
        cycle_warnings = []
        for smi in precursors:
            can = canonical(smi) or smi
            if can in self._seen:
                cycle_warnings.append(f"前体 {can} 已在树中，可能形成环路")

        # ── 写入树 ──
        evidence = TemplateEvidence(
            template_id=actual_template_id,
            template_name=template_name or fgi_template_name,
            reaction_category=actual_reaction_type,
            bond_atoms=list(bond) if bond else [],
            is_fgi=fgi_template_id is not None,
        )

        rxn = self.tree.add_reaction(
            product_smiles=smiles,
            precursors=precursors,
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

        for smi in precursors:
            can = canonical(smi) or smi
            is_term, cs_result = self._check_terminal(can)
            self.tree.update_complexity(can, cs_result)

            if is_term:
                self.tree.mark_terminal(can)
                new_terminal.append(can)
            else:
                self.tree.mark_intermediate(can)
                if can not in self._seen:
                    if len(self._queue) >= self.max_queue_size:
                        logger.warning("Queue full, marking %s as terminal", can)
                        self.tree.mark_terminal(can)
                        new_terminal.append(can)
                    else:
                        self._seen.add(can)
                        self._queue.append((can, depth + 1))
                        new_pending.append(can)
                else:
                    # 收敛路线: 已见过的分子不重复入队
                    logger.info("Convergent: %s already seen", can)

        self._current_context = None

        return CommitResult(
            success=True,
            reaction_node=rxn,
            new_pending=new_pending,
            new_terminal=new_terminal,
            tree_complete=self.is_complete() and self.tree.is_complete(),
            cycle_warnings=cycle_warnings,
            forward_validation=fv_result,
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
            outcome="terminal",
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
        result = {
            "bond_idx": bond_idx,
            "atoms": b.get("atoms", []),
            "bond_type": b.get("bond_type", ""),
            "in_ring": b.get("in_ring", False),
            "heuristic_score": b.get("heuristic_score", 0),
            "alternatives": b.get("alternatives", []),
            "hint": "用 try_disconnection(bond_idx, alt_idx) 沙盒试断。",
        }
        # 包含 smart_capping（如果 build_decision_context 已生成）
        if "smart_capping" in b:
            result["smart_capping"] = b["smart_capping"]
        else:
            # 按需生成
            try:
                from Rachel.chem_tools.smart_cap import suggest_capping
                cap_result = suggest_capping(ctx.smiles, tuple(b.get("atoms", [])))
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
            "hint": "用 try_fgi(fgi_idx) 沙盒试 FGI。",
        }

    # ── 沙盒试探: 不写入树 ──

    def _sandbox_classify_precursors(self, precursors: List[str]) -> Tuple[
        List[Dict[str, Any]], List[str]
    ]:
        """对前体做 CS 评分和 terminal 判定（不写入树）。"""
        details = []
        cycle_warnings = []
        for smi in precursors:
            can = canonical(smi) or smi
            cs = compute_cs_score(can)
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
                cycle_warnings.append(f"前体 {can} 已在树中")
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

        if bond is None:
            bonds = ctx.decision_context.get("disconnectable_bonds", [])
            if bond_idx < 0 or bond_idx >= len(bonds):
                return SandboxResult(
                    success=False,
                    error=f"bond_idx {bond_idx} out of range",
                )
            b = bonds[bond_idx]
            actual_bond = tuple(b["atoms"])
            alts = b.get("alternatives", [])
            if alt_idx < 0 or alt_idx >= len(alts):
                return SandboxResult(
                    success=False,
                    error=f"alt_idx {alt_idx} out of range (0-{len(alts)-1})",
                )
            alt = alts[alt_idx]
            actual_type = alt.get("template", "").split("(")[0].strip()
            actual_tid = alt.get("template_id", "")
            actual_tname = alt.get("template", "")

        # 执行断键
        break_result = execute_disconnection(smiles, actual_bond, actual_type)
        if not break_result.success:
            return SandboxResult(
                success=False,
                error=break_result.error or "disconnection failed",
                reaction_type=actual_type,
                template_name=actual_tname,
            )

        precursors = break_result.precursors
        if not precursors:
            return SandboxResult(success=False, error="no precursors generated")

        # 前体分析
        details, cycle_warnings = self._sandbox_classify_precursors(precursors)

        # 正向验证
        fv = None
        try:
            fv = validate_forward(
                precursors, smiles,
                template_id=actual_tid or None,
                reaction_category=actual_type or None,
            )
        except Exception:
            pass

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
            atom_balance=ab,
            cycle_warnings=cycle_warnings,
            reaction_type=actual_type,
            template_id=actual_tid,
            template_name=actual_tname,
        )

    def try_fgi(self, fgi_idx: int = 0) -> SandboxResult:
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

        fgi_result = execute_fgi(smiles, tid)
        if not fgi_result.success:
            return SandboxResult(
                success=False,
                error=fgi_result.error or "FGI failed",
                template_name=tname,
            )

        precursors = fgi_result.precursors
        details, cycle_warnings = self._sandbox_classify_precursors(precursors)

        fv = None
        try:
            fv = validate_forward(precursors, smiles, template_id=tid)
        except Exception:
            pass

        return SandboxResult(
            success=True,
            precursors=precursors,
            precursor_details=details,
            forward_validation=fv,
            cycle_warnings=cycle_warnings,
            reaction_type="fgi",
            template_id=tid,
            template_name=tname,
        )

    def try_precursors(
        self,
        precursors: List[str],
        reaction_type: str = "llm_proposed",
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
        valid_precursors = []
        for smi in parsed:
            can = canonical(smi)
            if not can:
                return SandboxResult(
                    success=False,
                    error=f"invalid SMILES: {smi}",
                )
            valid_precursors.append(can)

        # 前体分析
        details, cycle_warnings = self._sandbox_classify_precursors(valid_precursors)

        # 正向验证
        fv = None
        try:
            fv = validate_forward(
                valid_precursors, smiles,
                reaction_category=reaction_type,
            )
        except Exception:
            pass

        # 原子平衡
        ab = None
        try:
            ab = check_atom_balance(valid_precursors, smiles)
        except Exception:
            pass

        return SandboxResult(
            success=True,
            precursors=valid_precursors,
            precursor_details=details,
            forward_validation=fv,
            atom_balance=ab,
            cycle_warnings=cycle_warnings,
            reaction_type=reaction_type,
        )

    def finalize(self, llm_summary: str = "") -> Dict[str, Any]:
        """完成编排，返回最终 JSON。"""
        # 清空队列中剩余分子
        while self._queue:
            smi, _ = self._queue.popleft()
            self.tree.mark_terminal(smi)

        if self._current_context:
            self.tree.mark_terminal(self._current_context.smiles)
            self._current_context = None

        if self.tree.is_complete():
            self.tree.complete(llm_summary)
        else:
            self.tree.fail(llm_summary or "incomplete")

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
            print(f"完成: steps={status['steps_executed']}  "
                  f"molecules={status['total_molecules']}  "
                  f"depth={status['max_depth']}  "
                  f"elapsed={status['elapsed_sec']}s")

        return report
