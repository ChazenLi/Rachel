# Rachel/main — 多步逆合成编排引擎
# 合成树数据模型 + BFS 编排 + 审计状态 + 报告生成

from .retro_tree import (
    RetrosynthesisTree,
    MoleculeNode,
    ReactionNode,
    MoleculeRole,
    TreeStatus,
    LLMDecision,
    TemplateEvidence,
)

from .retro_state import (
    SynthesisAuditState,
    DecisionRecord,
    FailedAttempt,
    ProtectionEntry,
)

from .retro_orchestrator import (
    RetrosynthesisOrchestrator,
    ProposalContext,
    CommitResult,
    SandboxResult,
)

from .retro_report import (
    generate_forward_report,
    get_terminal_list,
    to_visualization_data,
)

from .retro_session import RetroSession
from .retro_cmd import RetroCmd
from .retro_output import export_results
from .retro_visualizer import generate_visual_report
