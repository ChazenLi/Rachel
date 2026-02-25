#!/usr/bin/env python
"""多步逆合成编排引擎 — 面向 LLM 的递归逆合成规划系统。

核心架构:
  - RetroNode: 合成树节点（分子 + 分析上下文 + 子节点）
  - RetroTree: 合成树管理（递归展开 + 回溯 + 状态追踪）
  - StepExecutor: 单步执行器（断键/FGI → 验证 → 评分）
  - MultiStepOrchestrator: 顶层编排器（LLM 交互循环 + 终止判断）

用法:
    # 1. 程序化调用（LLM 集成）
    from Rachel.tools.multistep_retro import MultiStepOrchestrator
    orch = MultiStepOrchestrator(max_depth=6, max_branches=3)
    tree = orch.init_tree("CC(=O)Oc1ccccc1C(=O)O")
    # LLM 循环: 获取待展开节点 → 生成上下文 → LLM 决策 → 执行
    pending = orch.get_pending_nodes()
    ctx = orch.get_node_context(pending[0])
    result = orch.execute_step(pending[0], bond=[3,5], reaction_type="ester_hydrolysis")
    # 重复直到所有叶节点 terminal

    # 2. 命令行演示（手动模拟 LLM 决策）
    python Rachel/tools/multistep_retro.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --auto

输出:
  - 完整合成树（JSON 可序列化）
  - 每步的验证结果和评分
  - 可视化就绪的树结构
  - 合成报告数据
"""

from __future__ import annotations

import json
import sys
import time
import uuid
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.chem_tools._rdkit_utils import canonical, parse_mol
from Rachel.chem_tools.mol_info import analyze_molecule
from Rachel.chem_tools.fg_detect import detect_functional_groups
from Rachel.chem_tools.cs_score import compute_cs_score, classify_complexity, score_progress
from Rachel.chem_tools.template_scan import find_disconnectable_bonds
from Rachel.chem_tools.bond_break import (
    execute_disconnection,
    execute_fgi,
    preview_disconnections,
)
from Rachel.chem_tools.forward_validate import validate_forward
from Rachel.chem_tools.fg_warnings import (
    check_fg_conflicts,
    suggest_protection_needs,
    check_reaction_compatibility,
)
from Rachel.tools.llm_retro_platform import build_decision_context, format_context_text


# ─────────────────────────────────────────────────────────────────────────
# 数据结构
# ─────────────────────────────────────────────────────────────────────────

class NodeStatus(str, Enum):
    """合成树节点状态。"""
    PENDING = "pending"          # 待分析（还没生成决策上下文）
    ANALYZED = "analyzed"        # 已分析（有决策上下文，等 LLM 决策）
    EXPANDED = "expanded"        # 已展开（LLM 已选断键，子节点已生成）
    TERMINAL = "terminal"        # 终端分子（可商购 / trivial）
    FAILED = "failed"            # 展开失败（无可行断键 / 验证不通过）


@dataclass
class StepResult:
    """单步逆合成执行结果。"""
    success: bool
    precursors: List[str] = field(default_factory=list)
    reaction_type: str = ""
    template_id: str = ""
    template_name: str = ""
    bond: Optional[Tuple[int, int]] = None
    is_fgi: bool = False
    validation: Optional[Dict[str, Any]] = None
    protection_notes: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    elapsed_sec: float = 0.0


@dataclass
class RetroNode:
    """合成树节点 — 代表一个分子及其逆合成状态。"""
    node_id: str
    smiles: str
    depth: int
    parent_id: Optional[str] = None
    status: NodeStatus = NodeStatus.PENDING

    # 分析数据（build_decision_context 的结果）
    decision_context: Optional[Dict[str, Any]] = None
    complexity: Optional[Dict[str, Any]] = None

    # 展开数据（LLM 决策 + 执行结果）
    step_result: Optional[StepResult] = None
    children_ids: List[str] = field(default_factory=list)

    # 元数据
    created_at: float = field(default_factory=time.time)
    analysis_time: float = 0.0
    expansion_time: float = 0.0

    @property
    def is_terminal(self) -> bool:
        return self.status == NodeStatus.TERMINAL

    @property
    def is_leaf(self) -> bool:
        return len(self.children_ids) == 0

    def to_dict(self) -> Dict[str, Any]:
        """序列化为字典（用于 JSON 输出和可视化）。"""
        d = {
            "node_id": self.node_id,
            "smiles": self.smiles,
            "depth": self.depth,
            "status": self.status.value,
            "parent_id": self.parent_id,
            "children_ids": self.children_ids,
        }
        if self.complexity:
            d["complexity"] = self.complexity
        if self.step_result:
            d["step"] = {
                "reaction_type": self.step_result.reaction_type,
                "template_name": self.step_result.template_name,
                "bond": list(self.step_result.bond) if self.step_result.bond else None,
                "is_fgi": self.step_result.is_fgi,
                "precursors": self.step_result.precursors,
                "validation_pass": (
                    self.step_result.validation.get("assessment", {}).get("pass")
                    if self.step_result.validation else None
                ),
                "feasibility_score": (
                    self.step_result.validation.get("assessment", {}).get("feasibility_score")
                    if self.step_result.validation else None
                ),
            }
        return d


# ─────────────────────────────────────────────────────────────────────────
# RetroTree — 合成树管理
# ─────────────────────────────────────────────────────────────────────────

class RetroTree:
    """管理整棵逆合成树的状态。"""

    def __init__(self, target_smiles: str):
        can = canonical(target_smiles)
        if can is None:
            raise ValueError(f"Invalid target SMILES: {target_smiles}")
        self.target_smiles = can
        self.nodes: Dict[str, RetroNode] = {}
        self.root_id: str = ""
        self._init_root()

    def _init_root(self):
        root = RetroNode(
            node_id=self._gen_id(),
            smiles=self.target_smiles,
            depth=0,
        )
        self.root_id = root.node_id
        self.nodes[root.node_id] = root

    @staticmethod
    def _gen_id() -> str:
        return uuid.uuid4().hex[:8]

    @property
    def root(self) -> RetroNode:
        return self.nodes[self.root_id]

    def add_child(self, parent_id: str, smiles: str) -> RetroNode:
        """为 parent 添加一个子节点（前体分子）。"""
        parent = self.nodes[parent_id]
        child = RetroNode(
            node_id=self._gen_id(),
            smiles=smiles,
            depth=parent.depth + 1,
            parent_id=parent_id,
        )
        self.nodes[child.node_id] = child
        parent.children_ids.append(child.node_id)
        return child

    def get_pending_nodes(self) -> List[RetroNode]:
        """获取所有待处理的叶节点（PENDING 或 ANALYZED 状态）。"""
        return [
            n for n in self.nodes.values()
            if n.status in (NodeStatus.PENDING, NodeStatus.ANALYZED)
        ]

    def get_leaves(self) -> List[RetroNode]:
        """获取所有叶节点。"""
        return [n for n in self.nodes.values() if n.is_leaf]

    def is_complete(self) -> bool:
        """检查合成树是否完成（所有叶节点都是 terminal 或 failed）。"""
        leaves = self.get_leaves()
        return all(
            n.status in (NodeStatus.TERMINAL, NodeStatus.FAILED)
            for n in leaves
        )

    def is_fully_resolved(self) -> bool:
        """检查合成树是否完全解析成功（所有叶节点都是 terminal）。"""
        leaves = self.get_leaves()
        return all(n.status == NodeStatus.TERMINAL for n in leaves)

    def get_depth(self) -> int:
        """当前树的最大深度。"""
        return max(n.depth for n in self.nodes.values()) if self.nodes else 0

    def get_path_to_root(self, node_id: str) -> List[RetroNode]:
        """从节点回溯到根的路径。"""
        path = []
        current = self.nodes.get(node_id)
        while current:
            path.append(current)
            current = self.nodes.get(current.parent_id) if current.parent_id else None
        return path

    def count_by_status(self) -> Dict[str, int]:
        """按状态统计节点数。"""
        counts: Dict[str, int] = {}
        for n in self.nodes.values():
            key = n.status.value
            counts[key] = counts.get(key, 0) + 1
        return counts

    def to_dict(self) -> Dict[str, Any]:
        """序列化整棵树。"""
        return {
            "target": self.target_smiles,
            "root_id": self.root_id,
            "total_nodes": len(self.nodes),
            "max_depth": self.get_depth(),
            "status_counts": self.count_by_status(),
            "complete": self.is_complete(),
            "fully_resolved": self.is_fully_resolved(),
            "nodes": {nid: n.to_dict() for nid, n in self.nodes.items()},
        }

    def to_nested_dict(self, node_id: Optional[str] = None) -> Dict[str, Any]:
        """序列化为嵌套树结构（用于可视化）。"""
        nid = node_id or self.root_id
        node = self.nodes[nid]
        d = node.to_dict()
        d["children"] = [
            self.to_nested_dict(cid) for cid in node.children_ids
        ]
        return d


# ─────────────────────────────────────────────────────────────────────────
# StepExecutor — 单步执行器
# ─────────────────────────────────────────────────────────────────────────

class StepExecutor:
    """执行单步逆合成：断键/FGI → 验证 → 评分。"""

    @staticmethod
    def analyze_node(
        node: RetroNode,
        target_smiles: str,
        terminal_cs_threshold: Optional[float] = None,
    ) -> None:
        """分析节点分子，生成决策上下文。就地修改 node。

        terminal 判定逻辑（三层）:
          1. 重原子数 ≤ 6 → 无条件 terminal（甲苯、苯胺级别）
          2. CS score ≤ terminal_cs_threshold → terminal
          3. 无可断键位 → terminal（工具层找不到拆法）
        """
        t0 = time.time()

        # 复杂度评分
        cs_result = compute_cs_score(node.smiles)
        node.complexity = cs_result

        # 层级 1: 极小分子无条件 terminal
        mol = parse_mol(node.smiles)
        if mol and mol.GetNumHeavyAtoms() <= 6:
            node.status = NodeStatus.TERMINAL
            node.analysis_time = time.time() - t0
            return

        # 层级 2: CS score 阈值判定
        threshold = terminal_cs_threshold if terminal_cs_threshold is not None else 2.5
        cs_score = cs_result.get("cs_score", 0)
        if cs_score <= threshold:
            node.status = NodeStatus.TERMINAL
            node.analysis_time = time.time() - t0
            return

        # 生成完整决策上下文
        ctx = build_decision_context(node.smiles)
        node.decision_context = ctx

        # 层级 3: 无可断键位 → terminal（没法再拆了）
        bonds = ctx.get("disconnectable_bonds", [])
        has_viable = any(b.get("alternatives") for b in bonds)
        if not bonds or not has_viable:
            node.status = NodeStatus.TERMINAL
            node.analysis_time = time.time() - t0
            return

        node.status = NodeStatus.ANALYZED
        node.analysis_time = time.time() - t0

    @staticmethod
    def execute_bond_disconnection(
        smiles: str,
        bond: Tuple[int, int],
        reaction_type: str,
        template_id: Optional[str] = None,
        target_smiles: Optional[str] = None,
    ) -> StepResult:
        """执行断键操作并验证。"""
        t0 = time.time()

        # 执行断键
        break_result = execute_disconnection(smiles, bond, reaction_type)

        if not break_result.success:
            return StepResult(
                success=False,
                bond=bond,
                reaction_type=reaction_type,
                error=break_result.error or "disconnection failed",
                elapsed_sec=time.time() - t0,
            )

        precursors = break_result.precursors

        # 正向验证
        validation = validate_forward(
            precursors, smiles,
            template_id=template_id,
            reaction_category=reaction_type,
        )

        # 保护基建议
        protection = None
        try:
            protection = suggest_protection_needs(smiles, reaction_type)
        except Exception:
            pass

        return StepResult(
            success=True,
            precursors=precursors,
            reaction_type=reaction_type,
            template_id=template_id or "",
            bond=bond,
            is_fgi=False,
            validation=validation,
            protection_notes=protection,
            elapsed_sec=time.time() - t0,
        )

    @staticmethod
    def execute_fgi_step(
        smiles: str,
        template_id: str,
        template_name: str = "",
    ) -> StepResult:
        """执行 FGI（官能团互变）操作并验证。"""
        t0 = time.time()

        fgi_result = execute_fgi(smiles, template_id)

        if not fgi_result.success:
            return StepResult(
                success=False,
                is_fgi=True,
                template_id=template_id,
                error=fgi_result.error or "FGI failed",
                elapsed_sec=time.time() - t0,
            )

        precursors = fgi_result.precursors

        # 正向验证
        validation = validate_forward(precursors, smiles, template_id=template_id)

        return StepResult(
            success=True,
            precursors=precursors,
            reaction_type="fgi",
            template_id=template_id,
            template_name=template_name,
            is_fgi=True,
            validation=validation,
            elapsed_sec=time.time() - t0,
        )

    @staticmethod
    def auto_select_best_bond(ctx: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """自动选择最佳断键方案（用于 auto 模式 / LLM fallback）。

        策略: 按 heuristic_score 降序，选第一个有 alternatives 的键位，
        取第一个 alternative。
        """
        bonds = ctx.get("disconnectable_bonds", [])
        if not bonds:
            return None

        # 按 heuristic_score 降序排序
        sorted_bonds = sorted(
            bonds, key=lambda b: b.get("heuristic_score", 0), reverse=True
        )

        for bond in sorted_bonds:
            alts = bond.get("alternatives", [])
            if alts:
                best_alt = alts[0]
                return {
                    "bond": tuple(bond["atoms"]),
                    "reaction_type": best_alt.get("template", "").split("(")[0].strip(),
                    "template_id": best_alt.get("template_id", ""),
                    "template_name": best_alt.get("template", ""),
                    "precursors_preview": best_alt.get("precursors", []),
                }

        # fallback: 用第一个键位，无模板
        first = sorted_bonds[0]
        return {
            "bond": tuple(first["atoms"]),
            "reaction_type": (first.get("template_names", ["unknown"])[0]
                              if first.get("template_names") else "unknown"),
            "template_id": "",
            "template_name": "",
            "precursors_preview": [],
        }


# ─────────────────────────────────────────────────────────────────────────
# MultiStepOrchestrator — 顶层编排器
# ─────────────────────────────────────────────────────────────────────────

class MultiStepOrchestrator:
    """多步逆合成编排器 — 管理 LLM 交互循环。

    LLM 集成流程:
        orch = MultiStepOrchestrator()
        tree = orch.init_tree(target_smiles)

        while not orch.is_done():
            # 1. 获取待展开节点
            pending = orch.get_pending_nodes()
            if not pending:
                break

            node = pending[0]

            # 2. 获取 LLM 决策上下文（已格式化的文本）
            ctx_text = orch.get_node_context_text(node.node_id)
            # → 发给 LLM，LLM 返回决策

            # 3. 执行 LLM 的决策
            result = orch.execute_step(
                node.node_id,
                bond=(3, 5),
                reaction_type="ester_hydrolysis",
                template_id="R_001",
            )

            # 4. 检查结果，决定是否回溯
            if not result.success:
                orch.mark_failed(node.node_id, result.error)

        # 5. 生成报告
        report = orch.generate_report()
    """

    def __init__(
        self,
        max_depth: int = 8,
        max_branches: int = 5,
        max_total_nodes: int = 100,
        validation_hard_gate: bool = True,
        terminal_cs_threshold: Optional[float] = None,
    ):
        """
        参数:
            max_depth: 最大递归深度
            max_branches: 每个节点最大分支数（预留）
            max_total_nodes: 合成树最大节点数
            validation_hard_gate: 是否启用正向验证硬门控
            terminal_cs_threshold: 自定义 terminal 判定阈值。
                None = 使用 cs_score 默认阈值 (2.5)。
                设为更小的值（如 1.5）可强制更深的拆解。
        """
        self.max_depth = max_depth
        self.max_branches = max_branches
        self.max_total_nodes = max_total_nodes
        self.validation_hard_gate = validation_hard_gate
        self.terminal_cs_threshold = terminal_cs_threshold

        self.tree: Optional[RetroTree] = None
        self.executor = StepExecutor()
        self.step_history: List[Dict[str, Any]] = []
        self._start_time: float = 0.0

    # ── 初始化 ──

    def init_tree(self, target_smiles: str) -> RetroTree:
        """初始化合成树并分析根节点。"""
        self.tree = RetroTree(target_smiles)
        self._start_time = time.time()

        # 分析根节点
        self.executor.analyze_node(
            self.tree.root, self.tree.target_smiles,
            terminal_cs_threshold=self.terminal_cs_threshold,
        )
        return self.tree

    # ── 节点查询 ──

    def get_pending_nodes(self) -> List[RetroNode]:
        """获取所有待 LLM 决策的节点（已分析但未展开）。"""
        if not self.tree:
            return []
        return [
            n for n in self.tree.nodes.values()
            if n.status == NodeStatus.ANALYZED
        ]

    def get_all_pending(self) -> List[RetroNode]:
        """获取所有需要处理的节点（PENDING + ANALYZED）。"""
        if not self.tree:
            return []
        return self.tree.get_pending_nodes()

    def get_node_context(self, node_id: str) -> Optional[Dict[str, Any]]:
        """获取节点的决策上下文（dict 格式，供程序化使用）。"""
        if not self.tree:
            return None
        node = self.tree.nodes.get(node_id)
        if not node:
            return None

        # 如果还没分析，先分析
        if node.status == NodeStatus.PENDING:
            self.executor.analyze_node(
                node, self.tree.target_smiles,
                terminal_cs_threshold=self.terminal_cs_threshold,
            )

        if node.status == NodeStatus.TERMINAL:
            return {"smiles": node.smiles, "is_terminal": True}

        return node.decision_context

    def get_node_context_text(self, node_id: str) -> str:
        """获取节点的决策上下文（格式化文本，供 LLM prompt 使用）。"""
        ctx = self.get_node_context(node_id)
        if not ctx:
            return ""
        if ctx.get("is_terminal"):
            return f"分子 {ctx['smiles']} 是终端分子（可商购），无需继续拆解。"

        node = self.tree.nodes[node_id]

        # 添加多步上下文信息
        header_lines = []
        header_lines.append(f"═══ 多步逆合成 — 第 {node.depth + 1} 步 ═══")
        header_lines.append(f"目标产物: {self.tree.target_smiles}")
        if node.depth > 0:
            path = self.tree.get_path_to_root(node_id)
            path_smiles = [n.smiles for n in reversed(path)]
            header_lines.append(f"路径: {' → '.join(path_smiles)}")
        header_lines.append("")

        body = format_context_text(ctx)
        return "\n".join(header_lines) + body

    # ── 执行 ──

    def execute_step(
        self,
        node_id: str,
        bond: Optional[Tuple[int, int]] = None,
        reaction_type: str = "",
        template_id: Optional[str] = None,
        fgi_template_id: Optional[str] = None,
        fgi_template_name: str = "",
    ) -> StepResult:
        """执行 LLM 的逆合成决策。

        参数:
            node_id: 要展开的节点 ID
            bond: 断键原子对 (i, j)，与 fgi_template_id 二选一
            reaction_type: 反应类型
            template_id: 断键模板 ID
            fgi_template_id: FGI 模板 ID（如果是 FGI 操作）
            fgi_template_name: FGI 模板名称
        """
        if not self.tree:
            return StepResult(success=False, error="tree not initialized")

        node = self.tree.nodes.get(node_id)
        if not node:
            return StepResult(success=False, error=f"node {node_id} not found")

        if node.status not in (NodeStatus.ANALYZED, NodeStatus.PENDING):
            return StepResult(
                success=False,
                error=f"node {node_id} status is {node.status.value}, cannot expand",
            )

        # 如果还没分析，先分析
        if node.status == NodeStatus.PENDING:
            self.executor.analyze_node(
                node, self.tree.target_smiles,
                terminal_cs_threshold=self.terminal_cs_threshold,
            )
            if node.status == NodeStatus.TERMINAL:
                return StepResult(success=True, precursors=[])

        t0 = time.time()

        # 执行断键或 FGI
        if fgi_template_id:
            result = self.executor.execute_fgi_step(
                node.smiles, fgi_template_id, fgi_template_name
            )
        elif bond is not None:
            result = self.executor.execute_bond_disconnection(
                node.smiles, bond, reaction_type,
                template_id=template_id,
                target_smiles=self.tree.target_smiles,
            )
        else:
            return StepResult(success=False, error="must specify bond or fgi_template_id")

        # 验证硬门控
        if result.success and self.validation_hard_gate:
            val = result.validation
            if val and val.get("assessment", {}).get("pass") is False:
                reasons = val.get("assessment", {}).get("hard_fail_reasons", [])
                result = StepResult(
                    success=False,
                    precursors=result.precursors,
                    bond=result.bond,
                    reaction_type=result.reaction_type,
                    validation=val,
                    error=f"validation hard gate: {', '.join(reasons)}",
                    elapsed_sec=time.time() - t0,
                )

        # 更新节点状态
        if result.success:
            node.step_result = result
            node.status = NodeStatus.EXPANDED
            node.expansion_time = time.time() - t0

            # 为每个前体创建子节点并分析
            for prec_smi in result.precursors:
                can_smi = canonical(prec_smi)
                if not can_smi:
                    continue

                # 检查深度限制
                if node.depth + 1 >= self.max_depth:
                    child = self.tree.add_child(node_id, can_smi)
                    child.status = NodeStatus.FAILED
                    child.complexity = {"note": "max depth reached"}
                    continue

                # 检查总节点数限制
                if len(self.tree.nodes) >= self.max_total_nodes:
                    child = self.tree.add_child(node_id, can_smi)
                    child.status = NodeStatus.FAILED
                    child.complexity = {"note": "max nodes reached"}
                    continue

                # 检查是否已经在树中出现过（避免循环）
                existing = [
                    n for n in self.tree.nodes.values()
                    if n.smiles == can_smi and n.node_id != node_id
                ]
                if existing and any(e.status == NodeStatus.TERMINAL for e in existing):
                    # 已知 terminal，直接标记
                    child = self.tree.add_child(node_id, can_smi)
                    child.status = NodeStatus.TERMINAL
                    child.complexity = existing[0].complexity
                    continue

                child = self.tree.add_child(node_id, can_smi)
                # 立即分析子节点（判断是否 terminal）
                self.executor.analyze_node(
                    child, self.tree.target_smiles,
                    terminal_cs_threshold=self.terminal_cs_threshold,
                )

        # 记录历史
        self.step_history.append({
            "step": len(self.step_history) + 1,
            "node_id": node_id,
            "smiles": node.smiles,
            "depth": node.depth,
            "success": result.success,
            "reaction_type": result.reaction_type,
            "precursors": result.precursors,
            "error": result.error,
            "elapsed_sec": result.elapsed_sec,
        })

        return result

    def mark_failed(self, node_id: str, reason: str = ""):
        """手动标记节点为失败。"""
        if not self.tree:
            return
        node = self.tree.nodes.get(node_id)
        if node:
            node.status = NodeStatus.FAILED
            node.complexity = {"note": reason or "manually marked failed"}

    # ── 自动模式 ──

    def auto_expand_node(self, node_id: str) -> StepResult:
        """自动选择最佳方案展开节点（不需要 LLM）。"""
        ctx = self.get_node_context(node_id)
        if not ctx or ctx.get("is_terminal"):
            return StepResult(success=True, precursors=[])

        # 先尝试断键
        selection = self.executor.auto_select_best_bond(ctx)
        if selection:
            result = self.execute_step(
                node_id,
                bond=selection["bond"],
                reaction_type=selection["reaction_type"],
                template_id=selection.get("template_id"),
            )
            if result.success:
                return result

        # 断键失败，尝试 FGI
        fgi_options = ctx.get("fgi_options", [])
        for fgi in fgi_options:
            result = self.execute_step(
                node_id,
                fgi_template_id=fgi["template_id"],
                fgi_template_name=fgi.get("template", ""),
            )
            if result.success:
                return result

        # 全部失败
        self.mark_failed(node_id, "no viable disconnection found")
        return StepResult(success=False, error="no viable disconnection found")

    def auto_run(self, verbose: bool = True) -> Dict[str, Any]:
        """全自动运行多步逆合成（贪心策略，用于演示/测试）。

        每步选 heuristic_score 最高的断键方案，递归展开所有非 terminal 前体。
        """
        if not self.tree:
            raise RuntimeError("call init_tree() first")

        iteration = 0
        max_iterations = self.max_total_nodes  # 安全上限

        while iteration < max_iterations:
            iteration += 1

            # 获取待展开节点
            pending = self.get_pending_nodes()
            if not pending:
                # 检查是否有 PENDING 状态的节点需要先分析
                unanalyzed = [
                    n for n in self.tree.nodes.values()
                    if n.status == NodeStatus.PENDING
                ]
                if unanalyzed:
                    for n in unanalyzed:
                        self.executor.analyze_node(
                            n, self.tree.target_smiles,
                            terminal_cs_threshold=self.terminal_cs_threshold,
                        )
                    pending = self.get_pending_nodes()

                if not pending:
                    break

            # 选择深度最浅的节点优先展开（BFS 策略）
            node = min(pending, key=lambda n: n.depth)

            if verbose:
                print(f"\n[Step {iteration}] depth={node.depth}  "
                      f"smiles={node.smiles[:60]}...")

            result = self.auto_expand_node(node.node_id)

            if verbose:
                if result.success:
                    prec_str = " + ".join(result.precursors)
                    print(f"  → {result.reaction_type}: {prec_str[:80]}")
                    # 显示子节点状态
                    for cid in self.tree.nodes[node.node_id].children_ids:
                        child = self.tree.nodes[cid]
                        status_icon = "✓" if child.is_terminal else "…"
                        print(f"    {status_icon} {child.smiles[:50]}  "
                              f"[{child.status.value}]")
                else:
                    print(f"  ✗ {result.error}")

            if self.tree.is_complete():
                break

        if verbose:
            print(f"\n{'='*60}")
            counts = self.tree.count_by_status()
            print(f"完成: {json.dumps(counts, ensure_ascii=False)}")
            print(f"总耗时: {time.time() - self._start_time:.1f}s")

        return self.generate_report()

    # ── 状态查询 ──

    def is_done(self) -> bool:
        """合成树是否已完成。"""
        return self.tree.is_complete() if self.tree else True

    def get_progress(self) -> Dict[str, Any]:
        """获取当前进度摘要。"""
        if not self.tree:
            return {"initialized": False}
        counts = self.tree.count_by_status()
        return {
            "initialized": True,
            "target": self.tree.target_smiles,
            "total_nodes": len(self.tree.nodes),
            "max_depth": self.tree.get_depth(),
            "status_counts": counts,
            "complete": self.tree.is_complete(),
            "fully_resolved": self.tree.is_fully_resolved(),
            "elapsed_sec": round(time.time() - self._start_time, 1),
        }


    # ── 报告生成 ──

    def generate_report(self) -> Dict[str, Any]:
        """生成合成报告（结构化数据，供可视化和文本渲染）。"""
        if not self.tree:
            return {"error": "tree not initialized"}

        leaves = self.tree.get_leaves()
        terminal_leaves = [n for n in leaves if n.status == NodeStatus.TERMINAL]
        failed_leaves = [n for n in leaves if n.status == NodeStatus.FAILED]

        # 收集所有反应步骤
        reactions = []
        for node in self.tree.nodes.values():
            if node.step_result and node.step_result.success:
                reactions.append({
                    "step_node": node.node_id,
                    "target": node.smiles,
                    "depth": node.depth,
                    "reaction_type": node.step_result.reaction_type,
                    "template_name": node.step_result.template_name,
                    "is_fgi": node.step_result.is_fgi,
                    "precursors": node.step_result.precursors,
                    "validation_pass": (
                        node.step_result.validation.get("assessment", {}).get("pass")
                        if node.step_result.validation else None
                    ),
                    "feasibility_score": (
                        node.step_result.validation.get("assessment", {}).get("feasibility_score")
                        if node.step_result.validation else None
                    ),
                })

        # 收集所有终端分子（起始原料）
        starting_materials = list({n.smiles for n in terminal_leaves})

        report = {
            "target": self.tree.target_smiles,
            "success": self.tree.is_fully_resolved(),
            "partial": self.tree.is_complete() and not self.tree.is_fully_resolved(),
            "summary": {
                "total_nodes": len(self.tree.nodes),
                "total_steps": len(reactions),
                "max_depth": self.tree.get_depth(),
                "n_starting_materials": len(starting_materials),
                "n_terminal": len(terminal_leaves),
                "n_failed": len(failed_leaves),
                "elapsed_sec": round(time.time() - self._start_time, 1),
            },
            "starting_materials": starting_materials,
            "reactions": reactions,
            "step_history": self.step_history,
            "tree": self.tree.to_dict(),
            "tree_nested": self.tree.to_nested_dict(),
        }

        return report

    def format_report_text(self, report: Optional[Dict[str, Any]] = None) -> str:
        """将报告格式化为人类可读文本。"""
        if report is None:
            report = self.generate_report()

        lines = []
        lines.append("╔══════════════════════════════════════════════════╗")
        lines.append("║          多步逆合成分析报告                      ║")
        lines.append("╚══════════════════════════════════════════════════╝")
        lines.append("")

        # 目标
        lines.append(f"目标分子: {report['target']}")
        s = report["summary"]
        status = "✓ 完全解析" if report["success"] else (
            "△ 部分解析" if report["partial"] else "… 进行中"
        )
        lines.append(f"状态:     {status}")
        lines.append(f"总节点:   {s['total_nodes']}  反应步数: {s['total_steps']}  "
                     f"最大深度: {s['max_depth']}")
        lines.append(f"起始原料: {s['n_starting_materials']} 种  "
                     f"终端: {s['n_terminal']}  失败: {s['n_failed']}")
        lines.append(f"耗时:     {s['elapsed_sec']}s")

        # 起始原料
        if report["starting_materials"]:
            lines.append("")
            lines.append("── 起始原料 ──")
            for sm in report["starting_materials"]:
                lines.append(f"  • {sm}")

        # 反应步骤
        if report["reactions"]:
            lines.append("")
            lines.append("── 反应步骤 ──")
            for i, rxn in enumerate(report["reactions"]):
                prec_str = " + ".join(rxn["precursors"])
                val_str = ""
                if rxn["feasibility_score"] is not None:
                    val_str = f"  [可行性={rxn['feasibility_score']:.3f}]"
                fgi_tag = " (FGI)" if rxn["is_fgi"] else ""
                lines.append(
                    f"  {i+1}. {rxn['target'][:50]}"
                    f"\n     → {rxn['reaction_type']}{fgi_tag}: {prec_str[:70]}"
                    f"{val_str}"
                )

        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────
# LLM Prompt 生成辅助
# ─────────────────────────────────────────────────────────────────────────

def build_multistep_prompt(
    orch: MultiStepOrchestrator,
    node_id: str,
    system_instructions: str = "",
) -> str:
    """为 LLM 生成完整的多步逆合成决策 prompt。

    包含:
      - 系统指令（可选）
      - 当前合成树进度
      - 当前节点的决策上下文
      - 期望的输出格式
    """
    parts = []

    if system_instructions:
        parts.append(system_instructions)
        parts.append("")

    # 进度信息
    progress = orch.get_progress()
    parts.append(f"当前进度: 已展开 {progress['total_nodes']} 个节点, "
                 f"最大深度 {progress['max_depth']}, "
                 f"状态 {json.dumps(progress['status_counts'], ensure_ascii=False)}")
    parts.append("")

    # 节点上下文
    ctx_text = orch.get_node_context_text(node_id)
    parts.append(ctx_text)
    parts.append("")

    # 输出格式要求
    parts.append("═══ 请选择一个断键方案 ═══")
    parts.append("请以 JSON 格式回复:")
    parts.append(json.dumps({
        "action": "disconnect | fgi | skip",
        "bond": [0, 0],
        "reaction_type": "反应类型名称",
        "template_id": "模板ID（如有）",
        "reasoning": "选择理由（简要）",
    }, indent=2, ensure_ascii=False))

    return "\n".join(parts)


def parse_llm_decision(response_text: str) -> Optional[Dict[str, Any]]:
    """解析 LLM 返回的决策 JSON。

    容错处理: 尝试从文本中提取 JSON 块。
    """
    import re

    # 尝试直接解析
    try:
        return json.loads(response_text.strip())
    except json.JSONDecodeError:
        pass

    # 尝试提取 ```json ... ``` 块
    match = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", response_text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group(1))
        except json.JSONDecodeError:
            pass

    # 尝试提取第一个 { ... }
    match = re.search(r"\{[^{}]*\}", response_text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group(0))
        except json.JSONDecodeError:
            pass

    return None


def apply_llm_decision(
    orch: MultiStepOrchestrator,
    node_id: str,
    decision: Dict[str, Any],
) -> StepResult:
    """将 LLM 的决策应用到编排器。"""
    action = decision.get("action", "disconnect")

    if action == "skip":
        orch.mark_failed(node_id, "LLM chose to skip")
        return StepResult(success=False, error="LLM chose to skip")

    if action == "fgi":
        return orch.execute_step(
            node_id,
            fgi_template_id=decision.get("template_id", ""),
            fgi_template_name=decision.get("template_name", ""),
        )

    # disconnect (default)
    bond_raw = decision.get("bond", [])
    if not bond_raw or len(bond_raw) < 2:
        return StepResult(success=False, error="invalid bond in LLM decision")

    return orch.execute_step(
        node_id,
        bond=tuple(bond_raw[:2]),
        reaction_type=decision.get("reaction_type", ""),
        template_id=decision.get("template_id"),
    )


# ─────────────────────────────────────────────────────────────────────────
# CLI 入口
# ─────────────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(description="多步逆合成编排引擎")
    parser.add_argument("--smiles", type=str, required=True,
                        help="目标分子 SMILES")
    parser.add_argument("--auto", action="store_true",
                        help="全自动模式（贪心策略，不需要 LLM）")
    parser.add_argument("--max-depth", type=int, default=6,
                        help="最大递归深度 (default: 6)")
    parser.add_argument("--max-nodes", type=int, default=50,
                        help="最大节点数 (default: 50)")
    parser.add_argument("--terminal-cs", type=float, default=None,
                        help="Terminal CS 阈值覆盖。默认 None 使用三层判定。"
                             "设为 1.3 可强制多步拆解。")
    parser.add_argument("--json", action="store_true",
                        help="输出 JSON 格式报告")
    parser.add_argument("--save", type=str, default=None,
                        help="保存报告到文件")
    args = parser.parse_args()

    # 验证 SMILES
    mol = parse_mol(args.smiles)
    if mol is None:
        print(f"无效 SMILES: {args.smiles}")
        sys.exit(1)

    orch = MultiStepOrchestrator(
        max_depth=args.max_depth,
        max_total_nodes=args.max_nodes,
        terminal_cs_threshold=args.terminal_cs,
    )
    orch.init_tree(args.smiles)

    if args.auto:
        report = orch.auto_run(verbose=not args.json)
    else:
        # 交互模式: 输出第一步上下文，等待用户输入
        pending = orch.get_pending_nodes()
        if not pending:
            print("目标分子已是终端分子，无需拆解。")
            report = orch.generate_report()
        else:
            node = pending[0]
            print(orch.get_node_context_text(node.node_id))
            print("\n提示: 使用 --auto 进行全自动分析")
            report = orch.generate_report()

    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    elif args.auto:
        print("\n" + orch.format_report_text(report))

    if args.save:
        save_path = Path(args.save)
        save_path.write_text(
            json.dumps(report, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"\n报告已保存: {save_path}")


if __name__ == "__main__":
    main()
