# Rachel/chem_tools — 化学工具层

## 概述

`chem_tools` 是 Rachel 逆合成分析系统的化学计算工具层，为上层 LLM 逆合成决策引擎提供纯计算与事实输出支撑。本层严格遵循"只做计算和事实输出，不做决策提示"的核心原则——不生成任何 LLM prompt、不构建 context、不输出合成建议文本，仅返回结构化的化学计算结果。

本工具层由旧系统 `old-Rachel/chem_analyze`（25 个 Python 文件、60+ 公开函数、8 个类）重构迁移而来，精简为 **8 个 Python 模块、26 个公开函数、3 个 dataclass**，删除了 NetworkX 图层、LLM prompt 生成、compact 输出、druglikeness 判定、图像渲染等冗余代码，同时覆盖旧系统全部核心功能并新增智能断键推理（M8）。

---

## 目录结构

```
Rachel/chem_tools/
├── __init__.py              # 公开 API 导出（26 个函数）
├── _rdkit_utils.py          # M0: RDKit 公共工具（内部模块）
├── mol_info.py              # M1: 分子信息与图结构
├── fg_detect.py             # M2: 官能团识别与匹配
├── template_scan.py         # M3: 反应模板扫描与匹配
├── bond_break.py            # M4: 断键执行与前体生成（6 个函数）
├── forward_validate.py      # M5: 正向验证与原子守恒
├── cs_score.py              # M6: CS 合成难度评分
├── fg_warnings.py           # M7: 官能团排斥与保护警告
├── smart_cap.py             # M8: 智能断键推理（规则推断 capping 方案）
└── templates/               # 12 个 JSON 模板数据文件
    ├── reactions.json              # 369 retro + 311 forward 反应模板（680 条）
    ├── functional_groups.json      # 官能团 SMARTS 定义
    ├── fg_compatibility_matrix.json # 官能团-反应条件兼容性矩阵
    ├── protecting_groups.json      # 保护基定义、脱除条件、正交性
    ├── selectivity_conflicts.json  # 选择性冲突对
    ├── byproducts.json             # 副产物数据
    ├── reaction_role_requirements.json # 反应角色需求
    ├── reactive_sites.json         # 反应活性位点 SMARTS
    ├── structural_alerts.json      # 结构警报
    ├── dangerous_combos.json       # 危险官能团组合
    ├── known_scaffolds.json        # 已知骨架库
    └── selectivity_reactivity.json # 选择性-反应性分类
```

---

## 模块依赖关系

模块间维持严格的单向依赖（DAG），不存在循环依赖：

```
M0 (_rdkit_utils)          ← 被所有模块依赖
 ├── M1 (mol_info)
 ├── M2 (fg_detect)        ← 依赖 M0
 ├── M3 (template_scan)    ← 依赖 M0
 ├── M4 (bond_break)       ← 依赖 M0, M3
 ├── M5 (forward_validate) ← 依赖 M0, M2
 ├── M6 (cs_score)         ← 依赖 M0
 ├── M7 (fg_warnings)      ← 依赖 M0, M2
 └── M8 (smart_cap)        ← 依赖 M0
```

---

## 核心设计原则

1. 纯计算输出：所有函数只返回化学计算结果和事实数据，不包含任何 LLM prompt、决策建议或合成提示文本
2. 统一错误处理：所有公开函数对无效 SMILES 返回错误字典（如 `{"ok": False, "error": "invalid SMILES", "input": smiles}`），绝不抛出异常
3. 字典通信：模块间通过 `Dict[str, Any]` 进行数据通信，不使用自定义类作为模块间接口
4. 启发式标注：所有基于人工权重的评分输出均标注 `is_heuristic: True`，不伪装精确度
5. LRU 缓存加速：SMILES 解析（`parse_mol`）和 SMARTS 编译（`_compile_smarts`）均使用 `@lru_cache` 避免重复计算

---

## 公开 API 总览（26 个函数）

| 模块 | 函数 | 作用 |
|------|------|------|
| M1 | `analyze_molecule(smiles)` | 全面分子静态属性分析 |
| M1 | `match_known_scaffolds(smiles, max_hits)` | 已知骨架库匹配 |
| M2 | `detect_functional_groups(smiles)` | 官能团识别 |
| M2 | `detect_reactive_sites(smiles)` | 反应活性位点识别 |
| M2 | `detect_protecting_groups(smiles)` | 保护基检测 |
| M2 | `get_fg_reaction_mapping(groups)` | 官能团→反应类型映射 |
| M3 | `scan_applicable_reactions(smiles, mode)` | 反应模板扫描 |
| M3 | `find_disconnectable_bonds(smiles)` | 可断键位点识别 |
| M4 | `execute_disconnection(smiles, bond, reaction_type, custom_fragments)` | 单键断裂执行 |
| M4 | `execute_fgi(smiles, template_id)` | 官能团互变执行 |
| M4 | `execute_operations(smiles, operations, reaction_category)` | 操作序列执行 |
| M4 | `preview_disconnections(smiles)` | 预览所有可能的断键方案 |
| M4 | `try_retro_template(mol, bond, template_id)` | 指定模板逆合成执行 |
| M4 | `resolve_template_ids(reaction_name)` | 反应名称→模板 ID 解析 |
| M5 | `validate_forward(precursors, target, ...)` | 正向合成可行性验证 |
| M5 | `check_atom_balance(precursors, target, ...)` | 原子守恒检查 |
| M6 | `compute_cs_score(smiles)` | CS 合成难度评分 |
| M6 | `classify_complexity(smiles)` | 快速复杂度分类 |
| M6 | `score_progress(intermediate, target)` | 合成进度评估 |
| M6 | `batch_score_progress(intermediates, target)` | 批量进度评估 |
| M7 | `check_fg_conflicts(smiles)` | 官能团冲突检测 |
| M7 | `check_reaction_compatibility(smiles, reaction_category)` | 反应条件兼容性检查 |
| M7 | `suggest_protection_needs(smiles, planned_reaction)` | 保护基需求推荐 |
| M7 | `check_deprotection_safety(smiles, pg_to_remove)` | 脱保护安全性评估 |
| M8 | `suggest_capping(smiles, bond, ...)` | 智能断键 capping 推理 |
| M8 | `custom_cap(smiles, bond, cap_i, cap_j)` | 自定义 capping 基团 |

---

## 各模块详细说明

### M0: `_rdkit_utils.py` — RDKit 公共工具（内部模块）

集中管理所有 RDKit 基础操作，消除各模块重复的解析、匹配、相似度计算代码。本模块以下划线开头，不对外导出，仅供其他模块内部调用。

| 函数 | 原理 | 缓存 |
|------|------|------|
| `parse_mol(smiles)` | `Chem.MolFromSmiles` + `SanitizeMol` | `@lru_cache(1024)` |
| `canonical(smiles)` | `parse_mol` → `MolToSmiles` | 通过 parse_mol 缓存 |
| `validate_smiles(smiles)` | 5 步验证：非空→长度≤500→合法字符→嵌套深度≤20→RDKit 解析 | — |
| `smarts_match(mol, smarts)` | SMARTS 编译缓存 + `GetSubstructMatches` | `@lru_cache(512)` |
| `tanimoto(mol_a, mol_b)` | Morgan fingerprint (r=2, 2048-bit) → Tanimoto | — |
| `mol_formula_counter(mol)` | `AddHs` 后遍历原子统计 | — |
| `load_template(filename)` | JSON 加载 + 递归过滤 `__comment` | — |

常量：`ELECTRONEGATIVITY`（22 种元素）、`VALENCE_ELECTRONS`（22 种）、`FRAGMENT_TEMPLATES`（22 种常见片段）。

---

### M1: `mol_info.py` — 分子信息与图结构

`analyze_molecule(smiles)` 输出：smiles, formula, mw, descriptors, atoms, bonds, rings, stereo, symmetry, scaffold。

骨架拓扑标签：`fused_bicyclic`、`fused_polycyclic`、`bridged`、`spiro`、`macrocycle`（环≥12原子）、`strained_small_ring`（环≤4原子）、`simple_monocyclic`、`polycyclic_separated`、`acyclic`。

`match_known_scaffolds(smiles)` 三级匹配：SMARTS 子结构 → Murcko 精确 → Murcko 子结构。

---

### M2: `fg_detect.py` — 官能团识别与匹配

- `detect_functional_groups` — SMARTS 匹配，bridgehead 去重
- `detect_reactive_sites` — 按类别前缀分组（electrophilic/nucleophilic/...）
- `detect_protecting_groups` — 子结构包含去重 + 双向正交性检查
- `get_fg_reaction_mapping` — 70+ 条目映射表

---

### M3: `template_scan.py` — 反应模板扫描与匹配

`TemplateMatch` dataclass：template_id, name, category, type, prereq_smarts, rxn_smarts, matched_atoms, confidence, broken_bonds, incompatible_groups, conditions。

`find_disconnectable_bonds` 的 `heuristic_score = √(bond_priority × best_template_confidence)`，bond_priority: SINGLE=1.0, DOUBLE=0.8, TRIPLE=0.6, AROMATIC=0.4，环内键 ×0.7。

---

### M4: `bond_break.py` — 断键执行与前体生成

`BreakResult` dataclass：success, precursors, operation_log, error。

`execute_disconnection` 三级优先级回退：
1. custom_fragments（LLM 显式指定）
2. rxn SMARTS 模板执行
3. 简单断键 + 加 H（兜底）

三级 SanitizeMol 回退：标准 → 跳过属性检查 → 仅 Kekulize + SetAromaticity。

`execute_operations` 支持 6 种操作：disconnect, add_fragment, change_bond_order, ring_open, modify_atom, fused_ring_retro。

`execute_fgi` — 官能团互变执行。

`preview_disconnections` — 预览所有可能的断键方案及其前体。

---

### M5: `forward_validate.py` — 正向验证与原子守恒

`validate_forward` 5 步验证：

| 步骤 | 权重 | Hard Gate |
|------|------|-----------|
| 原子守恒 | 0.25 | severe/skeleton_imbalance |
| 模板正向执行 | 0.25 | — |
| MCS 骨架对齐 | 0.20 | aligned=False |
| 键变化拓扑 | 0.15 | has_hard_fail |
| 官能团兼容性 | 0.15 | forbidden FG |

`check_atom_balance` 采用两层小分子损失设计：

- Tier 1（类别特定）：根据 `reaction_category` 查 `_SMALL_MOL_LOSSES` 表（19 类反应），只允许该类反应已知的损失分子。
- Tier 2（通用回退）：无条件尝试 6 种非骨架小分子（H₂O, HCl, HBr, HI, HF, H₂）。不会自动解释 C/N/S 差异。

损失分子库 `_LOSS_COUNTERS` 共 13 种：H₂O, HCl, HBr, HI, HF, H₂, N₂, CO₂, MeOH, EtOH, AcOH, NH₃, SO₂。

不平衡检测为单向设计：
- `skeleton_imbalance`：仅当产物比前体多出 C/N/S 时触发（前体侧多出是正常的试剂离去基团，如 Wittig 的 Ph₃P=O）。
- `severe_imbalance`：仅当产物侧非 H 原子多出 > 4 时触发。

`balance_score` 为连续评分（0-1），公式 `1.0 - (unexplained_atoms / total_heavy_atoms)`，取代旧的二值 balanced/not 判定。直接参与 `feasibility_score` 加权平均（权重 0.25）。若触发 hard gate 则强制为 0.0。

返回字段包括：`balanced`, `balance_score`, `precursor_atoms`, `product_atoms`, `excess`, `deficit`, `adjusted_excess`, `adjusted_deficit`, `severe_imbalance`, `skeleton_imbalance`, `note`。

---

### M6: `cs_score.py` — CS 合成难度评分

6 维加权：size(0.55), ring(0.65), stereo(0.55), hetero(0.40), symmetry(-0.20), fg_density(0.35)。

范围 [1.0, 10.0]。分类：≤2.5 trivial, ≤6.0 moderate, >6.0 complex。

800 样本 USPTO50K 验证：中位数 4.21, p90 5.48。

---

### M7: `fg_warnings.py` — 官能团排斥与保护警告

- `check_fg_conflicts` — 选择性冲突 + 竞争官能团 + 危险组合
- `check_reaction_compatibility` — 反应家族兼容性矩阵
- `suggest_protection_needs` — 保护基推荐 + 正交性检查
- `check_deprotection_safety` — 脱保护条件安全性评估

---

### M8: `smart_cap.py` — 智能断键推理

`CappingProposal` dataclass：reaction_type, cap_i, cap_j, fragments, confidence, description。

`suggest_capping(smiles, bond)` 基于键两端原子环境（元素、杂化、芳香性、邻近官能团）推断 capping 方案，覆盖 13 类反应规则。

`custom_cap(smiles, bond, cap_i, cap_j)` 允许 LLM 自定义两端 capping 基团。

---

## 模板数据文件

| 文件 | 内容 | 使用模块 |
|------|------|---------|
| `reactions.json` | 680 条反应模板（369 retro + 311 forward） | M3, M4, M5 |
| `functional_groups.json` | 官能团 SMARTS 定义 | M2, M5, M7 |
| `fg_compatibility_matrix.json` | 反应家族→官能团兼容性矩阵 | M5, M7 |
| `protecting_groups.json` | 保护基定义 + 正交性映射 | M2, M7 |
| `selectivity_conflicts.json` | 选择性冲突对 | M7 |
| `selectivity_reactivity.json` | 官能团→选择性类别映射 | M7 |
| `dangerous_combos.json` | 危险官能团组合 | M7 |
| `reactive_sites.json` | 反应活性位点 SMARTS | M2 |
| `known_scaffolds.json` | 已知骨架库 | M1 |
| `byproducts.json` | 副产物数据 | 预留 |
| `reaction_role_requirements.json` | 反应角色需求 | 预留 |
| `structural_alerts.json` | 结构警报 | 预留 |

---

## 依赖

- Python ≥ 3.10
- RDKit（`rdkit-pypi` 或 conda 安装）
- 标准库：json, re, math, collections, functools, pathlib, dataclasses, logging, typing
