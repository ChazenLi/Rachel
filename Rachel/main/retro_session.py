"""
JSON 会话持久化层
=================
用一个 JSON 文件维护编排器的完整状态，解决三个核心问题:

1. SMILES 安全 — JSON 中 [N+](=O)[O-] 是纯字符串，不被转义
2. 上下文恢复 — LLM 对话断了，读 JSON 即可恢复
3. 沙盒可擦除 — sandbox 区域随时擦除重写，commit 后才进 tree

JSON 结构:
{
  "session_id": "abc123",
  "target": { smiles, name, cs_score },
  "config": { max_depth, max_steps, terminal_cs_threshold, ... },
  "status": { phase, steps_executed, pending_count, ... },

  "current": {                    ← LLM 当前工作区
    "smiles": "...",
    "node_id": "mol_3",
    "depth": 2,
    "analysis": { cs_score, classification, functional_groups, ... },
    "bond_summary": [...],        ← 精简键位概览
    "fgi_summary": [...],
    "sandbox": {                  ← 沙盒试探结果（可擦除）
      "attempts": [
        { precursors, validation, reaction_type, ... },
      ],
      "selected": null | index    ← LLM 选中的方案
    }
  },

  "tree": { ... },               ← RetrosynthesisTree.to_dict()
  "audit_state": { ... },        ← SynthesisAuditState.to_dict()
  "queue": [ [smiles, depth], ... ]
}

用法:
  session = RetroSession.create("CCO", name="Ethanol")
  session.save("path/to/session.json")

  # LLM 对话断了，恢复:
  session = RetroSession.load("path/to/session.json")
  orch = session.get_orchestrator()
"""

from __future__ import annotations

import json
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Rachel.chem_tools._rdkit_utils import canonical
from Rachel.tools.llm_retro_platform import build_decision_context
from .retro_orchestrator import (
    RetrosynthesisOrchestrator,
    ProposalContext,
    SandboxResult,
    CommitResult,
)
from .retro_tree import RetrosynthesisTree, MoleculeRole, LLMDecision
from .retro_state import SynthesisAuditState


class RetroSession:
    """JSON 持久化会话。包装 RetrosynthesisOrchestrator，
    每次操作后自动保存状态到 JSON 文件。

    核心设计:
      - 所有 SMILES 只存在 JSON 中，不经过终端/shell
      - sandbox 区域可随时擦除，commit 后才写入 tree
      - LLM 读 JSON 获取上下文，写 JSON 提交决策
      - 对话断了读 JSON 即可恢复完整状态
    """

    def __init__(
        self,
        orchestrator: RetrosynthesisOrchestrator,
        session_path: str,
        session_id: str = "",
    ):
        self.orch = orchestrator
        self.path = Path(session_path)
        self.session_id = session_id or uuid.uuid4().hex[:8]
        self._sandbox_attempts: List[Dict[str, Any]] = []
        self._sandbox_selected: Optional[int] = None

    # ── 工厂方法 ──

    @classmethod
    def create(
        cls,
        target_smiles: str,
        session_path: str = ".rachel_session.json",
        target_name: str = "",
        *,
        max_depth: int = 15,
        max_steps: int = 50,
        terminal_cs_threshold: float = 2.5,
    ) -> RetroSession:
        """创建新会话。"""
        orch = RetrosynthesisOrchestrator(
            target_smiles,
            target_name=target_name,
            max_depth=max_depth,
            max_steps=max_steps,
            terminal_cs_threshold=terminal_cs_threshold,
        )
        session = cls(orch, session_path)
        session.save()
        return session

    @classmethod
    def load(cls, session_path: str) -> RetroSession:
        """从 JSON 文件恢复会话。"""
        path = Path(session_path)
        if not path.exists():
            raise FileNotFoundError(f"Session file not found: {path}")

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        # 重建 orchestrator
        config = data.get("config", {})
        target = data.get("target", {})

        orch = RetrosynthesisOrchestrator(
            target["smiles"],
            target_name=target.get("name", ""),
            max_depth=config.get("max_depth", 15),
            max_steps=config.get("max_steps", 50),
            terminal_cs_threshold=config.get("terminal_cs_threshold", 2.5),
        )

        # 恢复树
        if "tree" in data:
            orch.tree = RetrosynthesisTree.from_dict(data["tree"])

        # 恢复审计状态
        if "audit_state" in data:
            orch.audit_state = SynthesisAuditState.from_dict(data["audit_state"])

        # 恢复队列
        from collections import deque
        orch._queue = deque()
        for item in data.get("queue", []):
            orch._queue.append(tuple(item))

        # 恢复计数器
        status = data.get("status", {})
        orch._steps_executed = status.get("steps_executed", 0)

        # 恢复 seen 集合
        orch._seen = set(data.get("seen_smiles", []))

        session = cls(
            orch, session_path,
            session_id=data.get("session_id", ""),
        )

        # 恢复 sandbox
        current = data.get("current") or {}
        sandbox = current.get("sandbox") or {}
        session._sandbox_attempts = sandbox.get("attempts", [])
        session._sandbox_selected = sandbox.get("selected")

        # 恢复 _current_context（关键：让 explore/try_* 跨进程可用）
        if current and current.get("decision_tier") == "standard":
            cur_smiles = current.get("smiles", "")
            cur_depth = current.get("depth", 0)
            if cur_smiles:
                # 重新生成 decision_context（包含完整断键方案）
                decision_ctx = build_decision_context(cur_smiles)

                node = orch.tree.get_molecule_by_smiles(cur_smiles)
                is_target = (node.role == MoleculeRole.TARGET.value) if node else False

                ctx = ProposalContext(
                    smiles=cur_smiles,
                    node_id=current.get("node_id", node.node_id if node else ""),
                    depth=cur_depth,
                    cs_score=current.get("cs_score", 0),
                    classification=current.get("classification", ""),
                    is_terminal=current.get("is_terminal", False),
                    is_target=is_target,
                    depth_limited=cur_depth > config.get("max_depth", 15),
                    decision_context=decision_ctx,
                    seen_smiles=set(data.get("seen_smiles", [])),
                    steps_executed=status.get("steps_executed", 0),
                    steps_remaining=max(0, config.get("max_steps", 50) - status.get("steps_executed", 0)),
                    audit_state_summary=orch.audit_state.get_summary(),
                    failed_attempts_for_current=orch.audit_state.get_failures_for_molecule(cur_smiles),
                    decision_tier="standard",
                )
                orch._current_context = ctx

        return session

    # ── 序列化 ──

    def to_dict(self) -> Dict[str, Any]:
        """导出完整会话状态。"""
        orch = self.orch
        target_node = orch.tree.get_molecule_by_smiles(orch.tree.target)

        data: Dict[str, Any] = {
            "session_id": self.session_id,
            "target": {
                "smiles": orch.tree.target,
                "name": orch.tree.target_name,
                "cs_score": target_node.cs_score if target_node else 0,
            },
            "config": {
                "max_depth": orch.max_depth,
                "max_steps": orch.max_steps,
                "max_queue_size": orch.max_queue_size,
                "terminal_cs_threshold": orch.terminal_cs_threshold,
                "auto_forward_validate": orch.auto_forward_validate,
            },
            "status": orch.get_status(),
        }

        # 当前工作区
        ctx = orch._current_context
        if ctx and ctx.decision_tier != "quick_pass":
            current: Dict[str, Any] = {
                "smiles": ctx.smiles,
                "node_id": ctx.node_id,
                "depth": ctx.depth,
                "cs_score": ctx.cs_score,
                "classification": ctx.classification,
                "is_terminal": ctx.is_terminal,
                "is_target": ctx.is_target,
                "decision_tier": ctx.decision_tier,
            }

            # 精简键位概览
            if ctx.decision_context:
                dc = ctx.decision_context
                current["molecule"] = dc.get("molecule", {})
                current["functional_groups"] = dc.get("functional_groups", [])
                current["complexity"] = dc.get("complexity", {})
                current["warnings"] = dc.get("warnings", [])

                bonds = dc.get("disconnectable_bonds", [])
                current["bond_summary"] = [
                    {
                        "bond_idx": i,
                        "atoms": b.get("atoms", []),
                        "bond_type": b.get("bond_type", ""),
                        "in_ring": b.get("in_ring", False),
                        "heuristic_score": b.get("heuristic_score", 0),
                        "n_alternatives": len(b.get("alternatives", [])),
                        "reaction_types": [
                            a.get("template", "").split("(")[0].strip()
                            for a in b.get("alternatives", [])[:5]
                        ],
                    }
                    for i, b in enumerate(bonds)
                ]

                fgi = dc.get("fgi_options", [])
                if fgi:
                    current["fgi_summary"] = [
                        {"fgi_idx": i, "template": f.get("template", "")}
                        for i, f in enumerate(fgi)
                    ]

            # 沙盒
            current["sandbox"] = {
                "attempts": self._sandbox_attempts,
                "selected": self._sandbox_selected,
                "n_attempts": len(self._sandbox_attempts),
            }

            data["current"] = current
        elif ctx and ctx.decision_tier == "quick_pass":
            data["current"] = {
                "smiles": ctx.smiles,
                "cs_score": ctx.cs_score,
                "decision_tier": "quick_pass",
                "hint": "terminal 分子，直接 accept_terminal",
            }
        else:
            data["current"] = None

        # 队列
        data["queue"] = list(orch._queue)
        data["seen_smiles"] = sorted(orch._seen)

        # 树和审计状态
        data["tree"] = orch.tree.to_dict()
        data["audit_state"] = orch.audit_state.to_dict()

        return data

    def save(self) -> Path:
        """保存到 JSON 文件。"""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        data = self.to_dict()
        with open(self.path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        return self.path

    # ── LLM 上下文输出（分层） ──

    def get_context(self, detail: str = "compact") -> Dict[str, Any]:
        """获取 LLM 决策上下文。

        detail:
          "status"  — 只看状态（最省 token）
          "compact" — 当前分子 + 键位概览 + 沙盒历史
          "full"    — 包含完整断键方案（慎用）
          "tree"    — 包含完整树（调试用）
        """
        ctx = self.orch._current_context
        if ctx is None:
            return {
                "status": self.orch.get_status(),
                "current": None,
                "hint": "调用 prepare_next() 获取下一个分子",
            }

        if detail == "status":
            return {"status": self.orch.get_status()}

        result: Dict[str, Any] = {
            "status": self.orch.get_status(),
            "current": ctx.to_dict(detail=detail),
        }

        # 附加沙盒历史
        if self._sandbox_attempts:
            result["sandbox_history"] = {
                "n_attempts": len(self._sandbox_attempts),
                "attempts": self._sandbox_attempts,
                "selected": self._sandbox_selected,
            }

        if detail == "tree":
            result["tree_text"] = self.orch.tree.print_tree()

        return result

    # ── 编排操作（自动保存） ──

    def prepare_next(self) -> Optional[Dict[str, Any]]:
        """取下一个分子，返回精简上下文 dict。自动保存。"""
        self._sandbox_attempts = []
        self._sandbox_selected = None

        ctx = self.orch.prepare_next()
        if ctx is None:
            self.save()
            return None

        # quick_pass 自动处理
        if ctx.decision_tier == "quick_pass":
            result = {
                "action": "auto_terminal",
                "smiles": ctx.smiles,
                "cs_score": ctx.cs_score,
                "classification": ctx.classification,
            }
            self.orch.accept_terminal(reason="quick_pass")
            self.save()
            return result

        self.save()
        return self.get_context("compact")

    def explore_bond(self, bond_idx: int) -> Dict[str, Any]:
        """按需展开键位详情。不保存（只读操作）。"""
        return self.orch.explore_bond(bond_idx)

    def explore_fgi(self) -> Dict[str, Any]:
        """按需展开 FGI 详情。不保存。"""
        return self.orch.explore_fgi()

    def sandbox_try(
        self,
        *,
        bond_idx: int = -1,
        alt_idx: int = 0,
        fgi_idx: int = -1,
        precursors: Optional[List[str]] = None,
        reaction_type: str = "",
    ) -> Dict[str, Any]:
        """统一沙盒试探入口。结果追加到 sandbox.attempts，自动保存。

        三种模式:
          precursors=[...]           → LLM 自提前体
          bond_idx=0, alt_idx=0      → 试断键
          fgi_idx=0                  → 试 FGI
        """
        sb: SandboxResult

        if precursors is not None:
            sb = self.orch.try_precursors(precursors, reaction_type)
        elif bond_idx >= 0:
            sb = self.orch.try_disconnection(bond_idx, alt_idx)
        elif fgi_idx >= 0:
            sb = self.orch.try_fgi(fgi_idx)
        else:
            return {"error": "specify precursors, bond_idx, or fgi_idx"}

        attempt = sb.to_dict()
        attempt["attempt_idx"] = len(self._sandbox_attempts)
        if precursors is not None:
            attempt["source"] = "llm_proposed"
            attempt["input_precursors"] = precursors
        elif bond_idx >= 0:
            attempt["source"] = f"bond[{bond_idx}].alt[{alt_idx}]"
        else:
            attempt["source"] = f"fgi[{fgi_idx}]"

        # ── 自动保护基建议：forbidden_fg 时附带 actionable 建议 ──
        fv = attempt.get("forward_validation", {})
        hard_fails = fv.get("hard_fail_reasons") or []
        if "forbidden_fg" in hard_fails and sb.success:
            try:
                from Rachel.chem_tools.fg_warnings import suggest_protection_needs
                ctx = self.orch._current_context
                target_smi = ctx.smiles if ctx else ""
                rtype = attempt.get("reaction_type", "")
                if target_smi and rtype:
                    prot = suggest_protection_needs(target_smi, rtype)
                    if prot.get("needs_protection"):
                        attempt["protection_needed"] = prot["suggestions"]
                        attempt["hint"] = (
                            "forbidden_fg 检测到官能团冲突。建议先保护相关基团，"
                            "再重新 try_precursors。详见 protection_needed 字段。"
                        )
            except Exception:
                pass

        self._sandbox_attempts.append(attempt)
        self.save()
        return attempt

    def sandbox_clear(self) -> None:
        """擦除沙盒。"""
        self._sandbox_attempts = []
        self._sandbox_selected = None
        self.save()

    def sandbox_select(self, attempt_idx: int) -> Dict[str, Any]:
        """选中沙盒中的某个方案。"""
        if attempt_idx < 0 or attempt_idx >= len(self._sandbox_attempts):
            return {"error": f"attempt_idx {attempt_idx} out of range"}
        self._sandbox_selected = attempt_idx
        self.save()
        return self._sandbox_attempts[attempt_idx]

    def commit(
        self,
        attempt_idx: Optional[int] = None,
        *,
        reasoning: str = "",
        confidence: str = "medium",
        rejected: Optional[List[Dict[str, str]]] = None,
    ) -> Dict[str, Any]:
        """提交决策，写入树。自动保存。

        attempt_idx: 沙盒方案索引（默认用 sandbox_selected）
        """
        idx = attempt_idx if attempt_idx is not None else self._sandbox_selected
        if idx is None:
            return {"error": "no attempt selected. call sandbox_select(idx) first"}
        if idx < 0 or idx >= len(self._sandbox_attempts):
            return {"error": f"attempt_idx {idx} out of range"}

        attempt = self._sandbox_attempts[idx]
        if not attempt.get("success"):
            return {"error": "selected attempt was not successful"}

        precursors = attempt.get("precursors", [])
        reaction_type = attempt.get("reaction_type", "llm_proposed")
        template_id = attempt.get("template_id", "")
        template_name = attempt.get("template_name", "")

        decision = LLMDecision(
            selection_reasoning=reasoning,
            confidence=confidence,
            rejected_alternatives=rejected or [],
        )

        result = self.orch.commit_decision(
            precursor_override=precursors,
            reaction_type=reaction_type,
            template_id=template_id,
            template_name=template_name,
            llm_decision=decision,
        )

        # 清空沙盒
        self._sandbox_attempts = []
        self._sandbox_selected = None
        self.save()

        return result.to_dict()

    def accept_terminal(self, reason: str = "") -> None:
        """标记当前分子为 terminal。自动保存。"""
        self.orch.accept_terminal(reason=reason)
        self._sandbox_attempts = []
        self._sandbox_selected = None
        self.save()

    def skip(self, reason: str = "") -> None:
        """跳过当前分子。自动保存。"""
        self.orch.skip_current(reason)
        self._sandbox_attempts = []
        self._sandbox_selected = None
        self.save()

    def finalize(self, summary: str = "") -> Dict[str, Any]:
        """完成编排。自动保存。"""
        result = self.orch.finalize(summary)
        self.save()
        return result

    def get_orchestrator(self) -> RetrosynthesisOrchestrator:
        """获取底层编排器（高级用法）。"""
        return self.orch
