"""
合成树数据模型
==============
双类型节点图：MoleculeNode（分子）+ ReactionNode（反应）。
canonical SMILES 去重，完整 JSON 序列化。

从 old-Rachel/graph/retrosynthesis_tree.py 迁移重构，
适配新工具层 chem_tools 的输出格式。
"""

from __future__ import annotations

import json
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, List, Optional

from Rachel.chem_tools._rdkit_utils import canonical


# ─────────────────────────────────────────────────────────────────────────
# 枚举
# ─────────────────────────────────────────────────────────────────────────

class MoleculeRole(str, Enum):
    TARGET = "target"
    INTERMEDIATE = "intermediate"
    TERMINAL = "terminal"


class TreeStatus(str, Enum):
    IN_PROGRESS = "in_progress"
    COMPLETE = "complete"
    FAILED = "failed"


# ─────────────────────────────────────────────────────────────────────────
# 反应节点子结构
# ─────────────────────────────────────────────────────────────────────────

@dataclass
class TemplateEvidence:
    """工具层提供的模板证据（来自 find_disconnectable_bonds / preview_disconnections）。"""
    template_id: str = ""
    template_name: str = ""
    reaction_category: str = ""
    bond_atoms: List[int] = field(default_factory=list)
    heuristic_score: float = 0.0
    is_fgi: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> TemplateEvidence:
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class LLMDecision:
    """LLM 的选择理由（LLM 填充，编排器记录）。"""
    selection_reasoning: str = ""
    confidence: str = "medium"          # high / medium / low
    rejected_alternatives: List[Dict[str, str]] = field(default_factory=list)
    protection_needed: bool = False
    risk_assessment: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> LLMDecision:
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ─────────────────────────────────────────────────────────────────────────
# 核心节点
# ─────────────────────────────────────────────────────────────────────────

def _flatten_fv(fv: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten forward_validation dict for serialization.

    Handles both the original full dict from validate_forward()
    (with nested 'assessment' sub-dict) and the already-flattened
    dict from a previous to_dict() round-trip.
    """
    if "assessment" in fv:
        # Original full dict from validate_forward
        assess = fv["assessment"]
        return {
            "ok": fv.get("ok", False),
            "pass": assess.get("pass"),
            "feasibility_score": assess.get("feasibility_score"),
            "hard_fail_reasons": assess.get("hard_fail_reasons"),
        }
    else:
        # Already flattened (from a previous save/load round-trip)
        return {
            "ok": fv.get("ok", False),
            "pass": fv.get("pass"),
            "feasibility_score": fv.get("feasibility_score"),
            "hard_fail_reasons": fv.get("hard_fail_reasons"),
        }


@dataclass
class MoleculeNode:
    """分子节点。同一分子只出现一次（canonical SMILES 去重）。

    属性:
        smiles: canonical SMILES
        node_id: 唯一标识 (mol_0, mol_1, ...)
        role: target / intermediate / terminal
        depth: 在树中的深度（target=0）
        complexity: compute_cs_score 的结果
        decision_context: build_decision_context 的结果（仅非 terminal）
        llm_analysis: LLM 对分子的结构理解（可选）
    """
    smiles: str
    node_id: str = ""
    role: str = MoleculeRole.INTERMEDIATE.value
    depth: int = 0
    complexity: Optional[Dict[str, Any]] = None
    decision_context: Optional[Dict[str, Any]] = None
    llm_analysis: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        c = canonical(self.smiles)
        if c:
            self.smiles = c

    @property
    def is_terminal(self) -> bool:
        return self.role == MoleculeRole.TERMINAL.value

    @property
    def cs_score(self) -> float:
        if self.complexity:
            return self.complexity.get("cs_score", 0.0)
        return 0.0

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "smiles": self.smiles,
            "node_id": self.node_id,
            "role": self.role,
            "depth": self.depth,
        }
        if self.complexity:
            d["complexity"] = {
                "cs_score": self.complexity.get("cs_score", 0),
                "classification": self.complexity.get("classification", ""),
                "is_terminal": self.complexity.get("is_terminal", False),
            }
        if self.llm_analysis:
            d["llm_analysis"] = self.llm_analysis
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> MoleculeNode:
        return cls(
            smiles=d["smiles"],
            node_id=d.get("node_id", ""),
            role=d.get("role", MoleculeRole.INTERMEDIATE.value),
            depth=d.get("depth", 0),
            complexity=d.get("complexity"),
            llm_analysis=d.get("llm_analysis"),
        )


@dataclass
class ReactionNode:
    """反应节点。记录一步逆合成拆解的完整信息。"""
    step_id: str = ""
    depth: int = 0
    reaction_smiles: str = ""               # "A.B>>Product" 格式
    product_node: str = ""                  # 产物 MoleculeNode.node_id
    reactant_nodes: List[str] = field(default_factory=list)
    reaction_type: str = ""                 # 反应类型名称
    template_evidence: Optional[TemplateEvidence] = None
    llm_decision: Optional[LLMDecision] = None
    forward_validation: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "step_id": self.step_id,
            "depth": self.depth,
            "reaction_smiles": self.reaction_smiles,
            "product_node": self.product_node,
            "reactant_nodes": list(self.reactant_nodes),
            "reaction_type": self.reaction_type,
        }
        if self.template_evidence:
            d["template_evidence"] = self.template_evidence.to_dict()
        if self.llm_decision:
            d["llm_decision"] = self.llm_decision.to_dict()
        if self.forward_validation is not None:
            fv = self.forward_validation
            d["forward_validation"] = _flatten_fv(fv)
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> ReactionNode:
        te = d.get("template_evidence")
        ld = d.get("llm_decision")
        return cls(
            step_id=d.get("step_id", ""),
            depth=d.get("depth", 0),
            reaction_smiles=d.get("reaction_smiles", ""),
            product_node=d.get("product_node", ""),
            reactant_nodes=d.get("reactant_nodes", []),
            reaction_type=d.get("reaction_type", ""),
            template_evidence=TemplateEvidence.from_dict(te) if te else None,
            llm_decision=LLMDecision.from_dict(ld) if ld else None,
            forward_validation=d.get("forward_validation"),
        )


# ─────────────────────────────────────────────────────────────────────────
# 辅助函数
# ─────────────────────────────────────────────────────────────────────────

def build_reaction_smiles(precursors: List[str], product: str) -> str:
    """构造 'A.B>>Product' 格式的反应 SMILES。"""
    lhs = ".".join(precursors)
    return f"{lhs}>>{product}"


def parse_precursors(raw) -> List[str]:
    """统一解析前体格式。

    支持:
      - "A.B" → ["A", "B"]
      - ["A", "B"] → ["A", "B"]
      - ["A.B"] → ["A", "B"]
      - "A + B" → ["A", "B"]
    """
    if isinstance(raw, str):
        if "+" in raw:
            parts = [s.strip() for s in raw.split("+")]
        else:
            parts = [s.strip() for s in raw.split(".")]
        return [p for p in parts if p]

    if isinstance(raw, list):
        result = []
        for item in raw:
            if isinstance(item, str):
                result.extend(parse_precursors(item))
            else:
                result.append(str(item))
        return result

    return []


# ─────────────────────────────────────────────────────────────────────────
# 顶层容器
# ─────────────────────────────────────────────────────────────────────────

class RetrosynthesisTree:
    """逆合成树。双类型节点图，canonical SMILES 去重。

    核心不变量:
      - molecule_nodes 以 canonical SMILES 去重（_smiles_index 维护映射）
      - node_id 自增: mol_0, mol_1, ...
      - step_id 自增: rxn_1, rxn_2, ...
      - reaction_nodes 按添加顺序排列
    """

    def __init__(self, target_smiles: str, target_name: str = ""):
        can = canonical(target_smiles)
        if not can:
            raise ValueError(f"Invalid target SMILES: {target_smiles}")

        self.target: str = can
        self.target_name: str = target_name
        self.tree_id: str = uuid.uuid4().hex[:8]
        self.created_at: str = datetime.now(timezone.utc).isoformat()
        self.status: str = TreeStatus.IN_PROGRESS.value
        self.llm_summary: str = ""

        # 节点存储
        self.molecule_nodes: Dict[str, MoleculeNode] = {}
        self.reaction_nodes: List[ReactionNode] = []

        # 索引
        self._smiles_index: Dict[str, str] = {}  # canonical SMILES → node_id
        self._mol_counter: int = 0
        self._rxn_counter: int = 0

        # 创建 target 节点
        self._create_molecule(can, role=MoleculeRole.TARGET.value, depth=0)

    # ── 内部方法 ──

    def _next_mol_id(self) -> str:
        mid = f"mol_{self._mol_counter}"
        self._mol_counter += 1
        return mid

    def _next_rxn_id(self) -> str:
        self._rxn_counter += 1
        return f"rxn_{self._rxn_counter}"

    def _create_molecule(
        self, smiles: str, role: str = "intermediate", depth: int = 0,
    ) -> MoleculeNode:
        """创建或获取分子节点（canonical SMILES 去重）。"""
        can = canonical(smiles) or smiles
        if can in self._smiles_index:
            return self.molecule_nodes[self._smiles_index[can]]

        node = MoleculeNode(
            smiles=can,
            node_id=self._next_mol_id(),
            role=role,
            depth=depth,
        )
        self.molecule_nodes[node.node_id] = node
        self._smiles_index[can] = node.node_id
        return node

    # ── 查询 ──

    def get_molecule_by_smiles(self, smiles: str) -> Optional[MoleculeNode]:
        can = canonical(smiles) or smiles
        nid = self._smiles_index.get(can)
        return self.molecule_nodes.get(nid) if nid else None

    def get_molecule(self, node_id: str) -> Optional[MoleculeNode]:
        return self.molecule_nodes.get(node_id)

    def has_molecule(self, smiles: str) -> bool:
        can = canonical(smiles) or smiles
        return can in self._smiles_index

    # ── 修改 ──

    def mark_terminal(self, smiles: str) -> Optional[MoleculeNode]:
        node = self.get_molecule_by_smiles(smiles)
        if node:
            node.role = MoleculeRole.TERMINAL.value
        return node

    def mark_intermediate(self, smiles: str) -> Optional[MoleculeNode]:
        node = self.get_molecule_by_smiles(smiles)
        if node:
            node.role = MoleculeRole.INTERMEDIATE.value
        return node

    def update_complexity(self, smiles: str, cs_result: Dict[str, Any]) -> None:
        node = self.get_molecule_by_smiles(smiles)
        if node:
            node.complexity = cs_result

    def update_decision_context(self, smiles: str, ctx: Dict[str, Any]) -> None:
        node = self.get_molecule_by_smiles(smiles)
        if node:
            node.decision_context = ctx

    def update_llm_analysis(self, smiles: str, analysis: Dict[str, Any]) -> None:
        node = self.get_molecule_by_smiles(smiles)
        if node:
            node.llm_analysis = analysis

    def add_reaction(
        self,
        product_smiles: str,
        precursors: List[str],
        reaction_type: str = "",
        template_evidence: Optional[TemplateEvidence] = None,
        llm_decision: Optional[LLMDecision] = None,
        forward_validation: Optional[Dict[str, Any]] = None,
        depth: int = 0,
    ) -> ReactionNode:
        """添加一步反应，自动创建前体分子节点。"""
        product_node = self.get_molecule_by_smiles(product_smiles)
        if not product_node:
            product_node = self._create_molecule(product_smiles, depth=depth)

        reactant_ids = []
        for smi in precursors:
            can = canonical(smi) or smi
            child = self._create_molecule(can, depth=depth + 1)
            reactant_ids.append(child.node_id)

        rxn_smiles = build_reaction_smiles(precursors, product_smiles)

        rxn = ReactionNode(
            step_id=self._next_rxn_id(),
            depth=depth,
            reaction_smiles=rxn_smiles,
            product_node=product_node.node_id,
            reactant_nodes=reactant_ids,
            reaction_type=reaction_type,
            template_evidence=template_evidence,
            llm_decision=llm_decision,
            forward_validation=forward_validation,
        )
        self.reaction_nodes.append(rxn)
        return rxn

    # ── 状态查询 ──

    def get_terminal_smiles(self) -> List[str]:
        return [n.smiles for n in self.molecule_nodes.values()
                if n.role == MoleculeRole.TERMINAL.value]

    def get_pending_molecules(self) -> List[MoleculeNode]:
        """获取所有非 terminal、非 target 且未被拆解的中间体。"""
        expanded_products = {rxn.product_node for rxn in self.reaction_nodes}
        return [
            n for n in self.molecule_nodes.values()
            if n.role == MoleculeRole.INTERMEDIATE.value
            and n.node_id not in expanded_products
        ]

    @property
    def max_depth(self) -> int:
        if not self.molecule_nodes:
            return 0
        return max(n.depth for n in self.molecule_nodes.values())

    @property
    def total_steps(self) -> int:
        return len(self.reaction_nodes)

    def is_complete(self) -> bool:
        """所有叶节点都是 terminal（没有未拆解的中间体）。"""
        expanded_products = {rxn.product_node for rxn in self.reaction_nodes}
        for node in self.molecule_nodes.values():
            if node.role == MoleculeRole.TERMINAL.value:
                continue
            if node.node_id not in expanded_products:
                # target 也需要被拆解（除非它本身就是 terminal）
                return False
        return True

    def complete(self, llm_summary: str = "") -> None:
        self.status = TreeStatus.COMPLETE.value
        self.llm_summary = llm_summary

    def fail(self, reason: str = "") -> None:
        self.status = TreeStatus.FAILED.value
        self.llm_summary = reason

    # ── 序列化 ──

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target": self.target,
            "target_name": self.target_name,
            "tree_id": self.tree_id,
            "created_at": self.created_at,
            "status": self.status,
            "llm_summary": self.llm_summary,
            "total_steps": self.total_steps,
            "max_depth": self.max_depth,
            "molecule_nodes": {
                nid: n.to_dict() for nid, n in self.molecule_nodes.items()
            },
            "reaction_nodes": [r.to_dict() for r in self.reaction_nodes],
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=False)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> RetrosynthesisTree:
        tree = cls.__new__(cls)
        tree.target = data["target"]
        tree.target_name = data.get("target_name", "")
        tree.tree_id = data.get("tree_id", uuid.uuid4().hex[:8])
        tree.created_at = data.get("created_at", "")
        tree.status = data.get("status", TreeStatus.IN_PROGRESS.value)
        tree.llm_summary = data.get("llm_summary", "")

        tree.molecule_nodes = {}
        tree._smiles_index = {}
        tree._mol_counter = 0
        tree._rxn_counter = 0

        for nid, nd in data.get("molecule_nodes", {}).items():
            node = MoleculeNode.from_dict(nd)
            tree.molecule_nodes[nid] = node
            tree._smiles_index[node.smiles] = nid
            # 更新计数器
            try:
                num = int(nid.split("_")[1])
                tree._mol_counter = max(tree._mol_counter, num + 1)
            except (IndexError, ValueError):
                pass

        tree.reaction_nodes = [
            ReactionNode.from_dict(rd) for rd in data.get("reaction_nodes", [])
        ]
        for rxn in tree.reaction_nodes:
            try:
                num = int(rxn.step_id.split("_")[1])
                tree._rxn_counter = max(tree._rxn_counter, num)
            except (IndexError, ValueError):
                pass

        return tree

    @classmethod
    def from_json(cls, json_str: str) -> RetrosynthesisTree:
        return cls.from_dict(json.loads(json_str))

    # ── 文本渲染 ──

    def print_tree(self) -> str:
        """打印合成树的文本表示。"""
        lines: List[str] = []
        target_node = self.get_molecule_by_smiles(self.target)
        if target_node:
            self._print_subtree(target_node.node_id, lines, prefix="", is_last=True)
        return "\n".join(lines)

    def _print_subtree(
        self, node_id: str, lines: List[str],
        prefix: str = "", is_last: bool = True,
    ) -> None:
        node = self.molecule_nodes.get(node_id)
        if not node:
            return

        connector = "└── " if is_last else "├── "
        role_icon = {"target": "🎯", "terminal": "✓", "intermediate": "○"}
        icon = role_icon.get(node.role, "?")
        cs = f" CS={node.cs_score:.1f}" if node.complexity else ""
        lines.append(f"{prefix}{connector}{icon} {node.smiles[:50]}{cs} [{node.node_id}]")

        # 找到以此节点为产物的反应
        child_prefix = prefix + ("    " if is_last else "│   ")
        for rxn in self.reaction_nodes:
            if rxn.product_node == node_id:
                rxn_label = rxn.reaction_type or rxn.step_id
                lines.append(f"{child_prefix}  ↓ {rxn_label}")
                for i, rid in enumerate(rxn.reactant_nodes):
                    is_last_child = (i == len(rxn.reactant_nodes) - 1)
                    self._print_subtree(rid, lines, child_prefix, is_last_child)
