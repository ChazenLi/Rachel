# Rachel/chem_tools — 化学工具层公开 API
# 导出 M1-M7 全部 22 个公开函数

# M1: 分子信息（2个）
from .mol_info import analyze_molecule, match_known_scaffolds

# M2: 官能团识别（4个）
from .fg_detect import (
    detect_functional_groups,
    detect_reactive_sites,
    detect_protecting_groups,
    get_fg_reaction_mapping,
)

# M3: 模板扫描（2个）
from .template_scan import scan_applicable_reactions, find_disconnectable_bonds

# M4: 断键执行（6个）
from .bond_break import (
    execute_disconnection,
    execute_fgi,
    execute_operations,
    preview_disconnections,
    try_retro_template,
    resolve_template_ids,
)

# M5: 正向验证（2个）
from .forward_validate import validate_forward, check_atom_balance

# M6: CS 评分（4个）
from .cs_score import (
    compute_cs_score,
    classify_complexity,
    score_progress,
    batch_score_progress,
)

# M7: 官能团警告（4个）
from .fg_warnings import (
    check_fg_conflicts,
    check_reaction_compatibility,
    suggest_protection_needs,
    check_deprotection_safety,
)

# M8: 智能断键推理（2个）
from .smart_cap import suggest_capping, custom_cap

__all__ = [
    # M1
    "analyze_molecule",
    "match_known_scaffolds",
    # M2
    "detect_functional_groups",
    "detect_reactive_sites",
    "detect_protecting_groups",
    "get_fg_reaction_mapping",
    # M3
    "scan_applicable_reactions",
    "find_disconnectable_bonds",
    # M4
    "execute_disconnection",
    "execute_fgi",
    "execute_operations",
    "preview_disconnections",
    "try_retro_template",
    "resolve_template_ids",
    # M5
    "validate_forward",
    "check_atom_balance",
    # M6
    "compute_cs_score",
    "classify_complexity",
    "score_progress",
    "batch_score_progress",
    # M7
    "check_fg_conflicts",
    "check_reaction_compatibility",
    "suggest_protection_needs",
    "check_deprotection_safety",
    # M8
    "suggest_capping",
    "custom_cap",
]
