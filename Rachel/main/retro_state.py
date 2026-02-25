"""
合成审计状态
============
持久化 LLM 在逆合成规划过程中积累的结构化决策状态。

精简版，保留 old-Rachel 最有价值的部分:
  - DecisionRecord: 每步决策记录（自动追加）
  - FailedAttempt: 失败记录 + 经验教训
  - ProtectionEntry: 保护基生命周期追踪
  - strategic_plan: LLM 的全局战略规划

每次 prepare_next() 时注入摘要到返回 JSON，
LLM 即使上下文截断也能看到关键状态。
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class DecisionRecord:
    """单条决策历史记录。"""
    step_id: str = ""
    molecule: str = ""
    action: str = ""            # decide / propose / accept-terminal / skip
    reaction_name: str = ""
    reasoning_summary: str = ""
    outcome: str = ""           # committed / gate_failed / skipped / terminal
    confidence: str = "medium"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> DecisionRecord:
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class FailedAttempt:
    """单条失败记录。"""
    step_number: int = 0
    molecule: str = ""
    attempted_reaction: str = ""
    failure_reason: str = ""
    lesson: str = ""            # LLM 事后补充的经验教训

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> FailedAttempt:
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class ProtectionEntry:
    """保护基生命周期记录。"""
    functional_group: str = ""  # phenol, ketone, aldehyde ...
    position: str = ""          # C3-OH, ring A ketone ...
    protection: str = ""        # TBS ether, acetal, 延迟引入 ...
    install_step: str = ""      # rxn_2 或 planned
    remove_step: str = ""       # rxn_5 或 pending
    status: str = "planned"     # planned / installed / removed

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> ProtectionEntry:
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class SynthesisAuditState:
    """逆合成规划的全局审计状态。

    随编排器一起序列化，每次 prepare_next 时注入摘要给 LLM。
    """

    # 战略规划（LLM 初始化时填充）
    # 格式: {"key_disconnections": [...], "convergent_design": "...", ...}
    strategic_plan: Dict[str, Any] = field(default_factory=dict)

    # 保护基状态表
    protections: List[ProtectionEntry] = field(default_factory=list)

    # 决策历史
    decision_history: List[DecisionRecord] = field(default_factory=list)

    # 失败记录
    failed_attempts: List[FailedAttempt] = field(default_factory=list)

    # 线性步数监控
    linear_step_count: int = 0
    max_linear_target: int = 8
    target_cs_score: float = 0.0

    # ── 动态步数指导 ──

    def set_target_complexity(self, cs_score: float) -> None:
        """根据目标分子 CS score 动态设置建议最大线性步数。"""
        self.target_cs_score = cs_score
        if cs_score < 2.0:
            self.max_linear_target = 4
        elif cs_score < 3.0:
            self.max_linear_target = 8
        elif cs_score < 5.0:
            self.max_linear_target = 15
        else:
            self.max_linear_target = 25

    # ── 自动追加 ──

    def record_decision(
        self,
        step_id: str,
        molecule: str,
        action: str,
        reaction_name: str = "",
        reasoning_summary: str = "",
        outcome: str = "committed",
        confidence: str = "medium",
    ) -> None:
        self.decision_history.append(DecisionRecord(
            step_id=step_id,
            molecule=molecule,
            action=action,
            reaction_name=reaction_name,
            reasoning_summary=reasoning_summary[:120],
            outcome=outcome,
            confidence=confidence,
        ))

    def record_failure(
        self,
        molecule: str,
        attempted_reaction: str,
        failure_reason: str,
    ) -> None:
        self.failed_attempts.append(FailedAttempt(
            step_number=self.linear_step_count,
            molecule=molecule,
            attempted_reaction=attempted_reaction,
            failure_reason=failure_reason[:200],
        ))

    def get_failures_for_molecule(self, smiles: str) -> List[Dict[str, Any]]:
        return [f.to_dict() for f in self.failed_attempts if f.molecule == smiles]

    # ── LLM 主动更新 ──

    def update_strategic_plan(self, plan: Dict[str, Any]) -> None:
        self.strategic_plan = plan

    def upsert_protection(self, entry: Dict[str, Any]) -> None:
        """追加或更新保护基条目（按 functional_group + position 匹配）。"""
        fg = entry.get("functional_group", "")
        pos = entry.get("position", "")
        for existing in self.protections:
            if existing.functional_group == fg and existing.position == pos:
                for k, v in entry.items():
                    if k in ProtectionEntry.__dataclass_fields__ and v:
                        setattr(existing, k, v)
                return
        self.protections.append(ProtectionEntry.from_dict(entry))

    # ── 摘要输出 ──

    def get_summary(self, max_decisions: int = 8, max_failures: int = 5) -> Dict[str, Any]:
        """返回精简摘要，控制 token 量。"""
        summary: Dict[str, Any] = {}

        if self.strategic_plan:
            summary["strategic_plan"] = self.strategic_plan

        if self.protections:
            summary["protections"] = [p.to_dict() for p in self.protections]

        if self.decision_history:
            recent = self.decision_history[-max_decisions:]
            summary["recent_decisions"] = [d.to_dict() for d in recent]
            if len(self.decision_history) > max_decisions:
                summary["total_decisions"] = len(self.decision_history)

        if self.failed_attempts:
            recent = self.failed_attempts[-max_failures:]
            summary["recent_failures"] = [f.to_dict() for f in recent]
            if len(self.failed_attempts) > max_failures:
                summary["total_failures"] = len(self.failed_attempts)

        summary["linear_steps"] = self.linear_step_count
        summary["linear_target"] = self.max_linear_target
        if self.linear_step_count >= self.max_linear_target:
            summary["linear_warning"] = (
                f"线性步数 {self.linear_step_count} 已达建议上限 "
                f"{self.max_linear_target}，建议审查收敛替代方案"
            )

        return summary

    # ── 序列化 ──

    def to_dict(self) -> Dict[str, Any]:
        return {
            "strategic_plan": self.strategic_plan,
            "protections": [p.to_dict() for p in self.protections],
            "decision_history": [d.to_dict() for d in self.decision_history],
            "failed_attempts": [f.to_dict() for f in self.failed_attempts],
            "linear_step_count": self.linear_step_count,
            "max_linear_target": self.max_linear_target,
            "target_cs_score": self.target_cs_score,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> SynthesisAuditState:
        state = cls()
        state.strategic_plan = d.get("strategic_plan", {})
        state.protections = [
            ProtectionEntry.from_dict(p) for p in d.get("protections", [])
        ]
        state.decision_history = [
            DecisionRecord.from_dict(r) for r in d.get("decision_history", [])
        ]
        state.failed_attempts = [
            FailedAttempt.from_dict(f) for f in d.get("failed_attempts", [])
        ]
        state.linear_step_count = d.get("linear_step_count", 0)
        state.max_linear_target = d.get("max_linear_target", 8)
        state.target_cs_score = d.get("target_cs_score", 0.0)
        return state
