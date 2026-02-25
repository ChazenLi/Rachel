# Rachel 技术参考

## 1. RetroCmd 命令详情

### init

创建新会话，初始化合成树和编排器。

```python
cmd.execute("init", {
    "target": "CC(=O)Nc1ccc(O)cc1",  # 必填，目标分子 SMILES
    "name": "Paracetamol",            # 可选，名称
    "terminal_cs_threshold": 2.5,     # CS ≤ 此值判定为 terminal
    "max_steps": 50,                  # 最大反应步数
    "max_depth": 15,                  # 最大树深度
})
```

返回: `{ok, session_id, session_file, target, name}`

### next

取下一个待决策分子。自动跳过所有 quick_pass terminal（CS 低、重原子少、无可断键），
只在遇到 standard 分子时停下。

```python
ctx = cmd.execute("next")
# 返回 compact 上下文，或 {"action": "queue_empty"} 表示编排完成
```

返回的 compact 上下文结构:

```json
{
  "status": { "target", "status", "steps_executed", "pending_count", ... },
  "current": {
    "action": "awaiting_decision",
    "decision_tier": "standard",
    "smiles": "CC(=O)Nc1ccc(O)cc1",
    "node_id": "mol_0",
    "depth": 0,
    "cs_score": 1.3,
    "classification": "trivial",
    "is_terminal": true,
    "is_target": true,
    "functional_groups": [ { "name": "amide_generic", "count": 1, "atoms": [...] }, ... ],
    "bond_summary": [
      {
        "bond_idx": 0,
        "atoms": [3, 4],
        "bond_type": "SINGLE",
        "in_ring": false,
        "heuristic_score": 0.603,
        "n_alternatives": 3,
        "reaction_types": ["Chan Lam", "SNAr N-Nucleophile", "Ullmann N Arylation"],
        "smart_capping": [ ... ]
      }
    ],
    "fgi_summary": [ { "fgi_idx": 0, "template": "Swern (Retro, oxidation reduction)" } ],
    "audit_state_summary": { "linear_steps": 0, "linear_target": 4 }
  }
}
```

### context

获取当前分子的上下文，支持四级详细度。

```python
cmd.execute("context", {"detail": "compact"})  # compact / full / status / tree
```

### explore

展开某键位的完整前体方案（含 SMILES）。

```python
cmd.execute("explore", {"bond_idx": 1})
```

返回:

```json
{
  "bond_idx": 1,
  "atoms": [1, 3],
  "bond_type": "SINGLE",
  "heuristic_score": 0.603,
  "alternatives": [
    {
      "template": "Amide Bond Formation (Retro, ester amide)",
      "template_id": "G02_amide_retro",
      "precursors": ["CC(=O)Cl", "[NH]c1ccc(O)cc1"],
      "incompatible_with": [...]
    }
  ],
  "smart_capping": [
    {
      "reaction_type": "Amide bond formation",
      "fragments": ["CC(=O)O", "Nc1ccccc1"],
      "confidence": 0.9,
      "description": "RC(=O)-NR₂ → RC(=O)OH + HNR₂"
    }
  ]
}
```

### explore_fgi

展开 FGI（官能团互变）方案。

```python
cmd.execute("explore_fgi")
```

### try_bond

沙盒试断键。不写入树。

```python
cmd.execute("try_bond", {"bond_idx": 1, "alt_idx": 0})
```

### try_precursors

沙盒试 LLM 自提前体。验证 SMILES 合法性、CS 评分、正向验证、原子平衡、环路检测。

```python
cmd.execute("try_precursors", {
    "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
    "reaction_type": "Schotten-Baumann acylation"
})
```

返回:

```json
{
  "success": true,
  "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
  "precursor_details": [
    { "smiles": "CC(=O)Cl", "cs_score": 1.255, "classification": "trivial", "is_terminal": true }
  ],
  "forward_validation": {
    "ok": true, "pass": true, "feasibility_score": 0.705, "hard_fail_reasons": null
  },
  "atom_balance": { "balanced": false, "balance_score": 0.85, "note": "adjusted deficit: Cl=1" },
  "reaction_type": "Schotten-Baumann acylation",
  "attempt_idx": 0,
  "source": "llm_proposed"
}
```

### try_fgi

沙盒试官能团互变。

```python
cmd.execute("try_fgi", {"fgi_idx": 0})
```

### smart_cap

智能断键推理 — 基于键两端化学环境自动推断 capping 方案。

```python
# 按 bond_idx 查询（需要当前 context）
cmd.execute("smart_cap", {"bond_idx": 0})

# 直接指定 SMILES + 原子对
cmd.execute("smart_cap", {"smiles": "CC(=O)Nc1ccccc1", "bond": [1, 3]})

# 限制返回数量
cmd.execute("smart_cap", {"bond_idx": 0, "max": 3})
```

返回:

```json
{
  "ok": true,
  "proposals": [
    {
      "reaction_type": "Amide bond formation",
      "fragments": ["CC(=O)O", "Nc1ccccc1"],
      "confidence": 0.90,
      "description": "RC(=O)-NR₂ → RC(=O)OH + HNR₂"
    }
  ]
}
```

### custom_cap

LLM 自定义 capping — 指定断键两端各加什么基团。

```python
cmd.execute("custom_cap", {
    "bond_idx": 0,
    "cap_i": "Br",
    "cap_j": "B(O)O"
})
```

### sandbox_list

查看沙盒中所有方案的摘要。

```python
cmd.execute("sandbox_list")
# → { n_attempts, selected, attempts: [{idx, source, success, precursors, forward_pass}] }
```

### sandbox_clear

清空沙盒中所有方案。

```python
cmd.execute("sandbox_clear")
```

### select + commit

选中方案并提交到树。

```python
cmd.execute("select", {"idx": 1})
cmd.execute("commit", {
    "idx": 1,                          # 沙盒方案索引
    "reasoning": "酰氯法条件温和...",    # 决策理由
    "confidence": "high",               # high / medium / low
    "rejected": [                       # 被拒绝的方案
        {"method": "Ac2O route", "reason": "副产物处理"}
    ]
})
```

commit 返回:

```json
{
  "success": true,
  "step_id": "rxn_1",
  "reaction_smiles": "CC(=O)Cl.Nc1ccc(O)cc1>>CC(=O)Nc1ccc(O)cc1",
  "new_pending": [],
  "new_terminal": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
  "tree_complete": true,
  "forward_validation": { "ok": true, "pass": true, "feasibility_score": 0.705 }
}
```

### accept / skip

```python
cmd.execute("accept", {"reason": "简单商业试剂"})  # 标记为 terminal
cmd.execute("skip", {"reason": "无可行方案"})       # 跳过
```

### tree / status / finalize / report / export

```python
cmd.execute("tree")      # → {tree (文本), terminal_count, pending_count, terminals}
cmd.execute("status")    # → {target, status, steps_executed, pending_count, ...}
cmd.execute("finalize", {"summary": "..."})  # 完成编排
cmd.execute("report")    # → {report (文本), starting_materials}
cmd.execute("export", {"name": "Losartan", "output_dir": "..."})  # 导出结果
```

export 返回:

```json
{
  "output_dir": "output/20260220_235417_Losartan",
  "files": ["SYNTHESIS_REPORT.html", "SYNTHESIS_REPORT.md", "report.txt", ...],
  "n_files": 8,
  "n_images": 12,
  "html_report": "output/.../SYNTHESIS_REPORT.html",
  "visualization_ok": true,
  "summary": "Losartan: 5 步, 6 种起始原料, 可视化=✓"
}
```

---

## 2. 数据结构

### MoleculeNode

```python
@dataclass
class MoleculeNode:
    smiles: str              # canonical SMILES
    node_id: str             # "mol_0", "mol_1", ...
    role: str                # "target" / "intermediate" / "terminal" / "pending"
    depth: int               # 树中深度
    complexity: Dict         # {cs_score, classification, is_terminal}
    decision_context: Dict   # build_decision_context() 的完整输出
    llm_analysis: Dict       # LLM 的分析记录
```

### ReactionNode

```python
@dataclass
class ReactionNode:
    step_id: str             # "rxn_1", "rxn_2", ...
    depth: int
    reaction_smiles: str     # "A.B>>C"
    product_node: str        # mol_id
    reactant_nodes: List[str]  # [mol_id, ...]
    reaction_type: str
    template_evidence: TemplateEvidence
    llm_decision: LLMDecision
    forward_validation: Dict
```

### TemplateEvidence

```python
@dataclass
class TemplateEvidence:
    template_id: str = ""
    template_name: str = ""
    category: str = ""
    confidence: float = 0.0
    source: str = ""         # "template" / "llm_proposed" / "fgi"
```

### LLMDecision

```python
@dataclass
class LLMDecision:
    selection_reasoning: str       # 选择理由
    confidence: str                # "high" / "medium" / "low"
    rejected_alternatives: List[Dict]  # [{method, reason}]
    protection_needed: bool
    risk_assessment: str
```

### SynthesisAuditState

```python
@dataclass
class SynthesisAuditState:
    strategic_plan: Dict           # LLM 的全局战略
    protections: List[ProtectionEntry]  # 保护基生命周期
    decision_history: List[DecisionRecord]  # 每步决策记录
    failed_attempts: List[FailedAttempt]    # 失败记录
    linear_step_count: int
    max_linear_target: int         # 根据 CS score 动态设置
    target_cs_score: float
```

### DecisionRecord

```python
@dataclass
class DecisionRecord:
    step_id: str
    molecule: str
    action: str            # decide / propose / accept-terminal / skip
    reaction_name: str
    reasoning_summary: str
    outcome: str           # committed / gate_failed / skipped / terminal
    confidence: str        # high / medium / low
```

### ProtectionEntry

```python
@dataclass
class ProtectionEntry:
    functional_group: str  # phenol, ketone, aldehyde ...
    position: str          # C3-OH, ring A ketone ...
    protection: str        # TBS ether, acetal, 延迟引入 ...
    install_step: str      # rxn_2 或 planned
    remove_step: str       # rxn_5 或 pending
    status: str            # planned / installed / removed
```

---

## 3. 化学工具层 (chem_tools)

26 个函数，8 个模块。编排器内部调用，LLM 通常不需要直接使用。

| 模块 | 函数 | 作用 |
|------|------|------|
| M0 _rdkit_utils | `parse_mol`, `canonical`, `validate_smiles`, `smarts_match`, `tanimoto`, `mol_formula_counter`, `load_template` | RDKit 公共工具（内部） |
| M1 mol_info | `analyze_molecule`, `match_known_scaffolds` | 分子基础信息 |
| M2 fg_detect | `detect_functional_groups`, `detect_reactive_sites`, `detect_protecting_groups`, `get_fg_reaction_mapping` | 官能团识别 |
| M3 template_scan | `scan_applicable_reactions`, `find_disconnectable_bonds` | 模板扫描 |
| M4 bond_break | `execute_disconnection`, `execute_fgi`, `execute_operations`, `preview_disconnections`, `try_retro_template`, `resolve_template_ids` | 断键执行（6 个） |
| M5 forward_validate | `validate_forward`, `check_atom_balance` | 正向验证 |
| M6 cs_score | `compute_cs_score`, `classify_complexity`, `score_progress`, `batch_score_progress` | 复杂度评分 |
| M7 fg_warnings | `check_fg_conflicts`, `check_reaction_compatibility`, `suggest_protection_needs`, `check_deprotection_safety` | 官能团冲突 |
| M8 smart_cap | `suggest_capping`, `custom_cap` | 智能断键推理 |

### CS Score 分级

| 范围 | 分类 | 含义 |
|------|------|------|
| ≤ 2.5 | trivial | 简单试剂，可作为合成终点（is_terminal=true） |
| 2.5 - 6.0 | moderate | 中等复杂度，需要继续拆解 |
| > 6.0 | complex | 复杂分子，需要精心设计路线 |

### Smart Capping (M8)

`suggest_capping(smiles, bond)` 根据键两端原子的化学环境自动推断断键后的 capping 方案。

覆盖 13 类反应规则：

| 规则 | 键类型 | Cap (i端 / j端) | 置信度 |
|------|--------|-----------------|--------|
| Suzuki coupling | Ar-Ar | Br / B(OH)₂ | 0.92 |
| Negishi coupling | Ar-Ar | Br / ZnCl | 0.70 |
| Stille coupling | Ar-Ar | Br / SnMe₃ | 0.60 |
| Amide bond formation | C(=O)-N | OH / H | 0.90 |
| Amide (acid chloride) | C(=O)-N | Cl / H | 0.80 |
| Ester hydrolysis | C(=O)-O | OH / H | 0.88 |
| N-alkylation (SN2) | N-C(sp3) | H / Br | 0.82 |
| Reductive amination | N-C(sp3) | H / =O | 0.70 |
| Williamson ether | O-C(sp3) | H / Br | 0.78 |
| Buchwald-Hartwig | Ar-N | Br / H | 0.80 |
| SNAr/Ullmann | Ar-O | F / H | 0.65 |
| Heck | C(sp2)-C(sp3) | H / Br | 0.55 |
| Grignard | C-C (通用) | Br / MgBr | 0.45 |

集成位置：
- `build_decision_context` 中每个 bond 自动附带 `smart_capping` 字段
- `explore` 返回中包含 `smart_capping`
- compact 视图 `bond_summary` 中有概览
- `smart_cap` 命令可独立调用
- `custom_cap` 命令可 LLM 自定义 capping 基团

理论最大值约 8.75。800 样本 USPTO50K 验证：trivial 5.0%, moderate 91.6%, complex 3.4%；中位数 4.21。

基准分子参考：atorvastatin ≈ 4.6, cholesterol ≈ 6.2, morphine ≈ 6.8, taxol ≈ 8.2。

权重：size=0.55, ring=0.65, stereo=0.55, hetero=0.40, symmetry=-0.20, fg_density=0.35。size 维度基于重原子数（非 MW），fg_density 同时计数保护基和复杂度官能团（酰胺/磺酰胺/氨基甲酸酯/脲/酯等）。

### 正向验证 (forward_validation)

`validate_forward` 5 步验证，加权评分：

| 步骤 | 权重 | Hard Gate |
|------|------|-----------|
| 原子守恒 (check_atom_balance) | 0.25 | severe_imbalance / skeleton_imbalance |
| 模板正向执行 | 0.25 | — |
| MCS 骨架对齐 | 0.20 | aligned=False |
| 键变化拓扑 | 0.15 | has_hard_fail |
| 官能团兼容性 | 0.15 | forbidden FG |

`feasibility_score` 范围 0-1:
- \> 0.7: 高可行性
- 0.5 - 0.7: 中等，可接受
- < 0.5: 低可行性，需要审视

#### 原子守恒 — 两层小分子损失设计

`check_atom_balance` 采用两层（Tier）设计处理反应中的小分子损失：

- Tier 1（类别特定）：根据 `reaction_category` 查表，只允许该类反应已知的损失分子。覆盖 19 类反应（ester_hydrolysis, friedel_crafts, elimination, curtius 等）。
- Tier 2（通用回退）：无论是否有 category，都尝试 6 种非骨架小分子损失（H₂O, HCl, HBr, HI, HF, H₂）。约束：通用回退不会自动解释骨架原子（C/N/S）差异。

两层依次应用于 deficit（前体多出）和 excess（产物多出）两侧。

#### 单向不平衡检测

关键设计：`skeleton_imbalance` 和 `severe_imbalance` 均为单向检测，只在产物侧有多余原子时触发：

- `skeleton_imbalance`：产物比前体多出 C/N/S 原子（原子凭空出现）。前体侧多出 C/N/S 是正常的（Wittig 的 Ph₃P=O、Stille 的 Bu₃SnCl 等试剂离去基团）。
- `severe_imbalance`：产物侧非 H 原子多出 > 4。前体侧多出通过连续评分软惩罚。

#### 连续 balance_score

`balance_score` 取代了旧的二值 balanced/not 判定，计算公式：

```
balance_score = 1.0 - (unexplained_atoms / total_heavy_atoms)
```

其中 `unexplained_atoms` = 两层损失扣除后仍无法解释的原子总数。典型值：
- 完全守恒反应：1.0
- Wittig（Ph₃P=O 离去）：~0.77
- Suzuki（硼酸离去）：~0.81
- Stille（锡试剂离去）：~0.5
- 严重不平衡：→ 0.0（被 hard gate 拦截）

`balance_score` 直接作为原子守恒维度（权重 0.25）参与 `feasibility_score` 加权平均。若触发 severe_imbalance 或 skeleton_imbalance，该维度强制为 0.0。

#### check_atom_balance 返回字段

```json
{
  "balanced": true,
  "balance_score": 0.77,
  "precursor_atoms": {"C": 20, "H": 18, "O": 2, "P": 1},
  "product_atoms": {"C": 10, "H": 10, "O": 1},
  "excess": {},
  "deficit": {"C": 10, "H": 8, "O": 1, "P": 1},
  "adjusted_excess": {},
  "adjusted_deficit": {"C": 10, "O": 1, "P": 1},
  "severe_imbalance": false,
  "skeleton_imbalance": false,
  "note": "balanced after accounting for small-molecule loss"
}
```

#### hard_fail_reasons 可能的值

| 值 | 含义 |
|----|------|
| `severe_imbalance` | 产物侧非 H 原子多出 > 4 |
| `skeleton_imbalance` | 产物比前体多出骨架原子（C/N/S） |
| `scaffold_not_aligned` | MCS 骨架对齐不足 |
| `bond_topology_violation` | 键变化拓扑不可行 |
| `forbidden_fg` | 含有反应禁忌官能团 |

---

## 4. JSON 会话文件结构

```json
{
  "session_id": "abc123",
  "target": { "smiles": "...", "name": "...", "cs_score": 1.3 },
  "config": { "max_depth": 15, "max_steps": 50, "terminal_cs_threshold": 2.5 },
  "status": { "status": "in_progress", "steps_executed": 1, "pending_count": 0 },
  "current": {
    "smiles": "...",
    "node_id": "mol_0",
    "sandbox": {
      "attempts": [ { "precursors": [...], "success": true, ... } ],
      "selected": null,
      "n_attempts": 1
    }
  },
  "queue": [ ["SMILES", depth], ... ],
  "seen_smiles": ["..."],
  "tree": { "molecule_nodes": {...}, "reaction_nodes": [...] },
  "audit_state": {
    "strategic_plan": {...},
    "protections": [...],
    "decision_history": [...],
    "failed_attempts": [...],
    "linear_step_count": 0,
    "max_linear_target": 8,
    "target_cs_score": 0.0
  }
}
```

---

## 5. 终止判定逻辑

三层判定（任一满足即为 terminal）:

1. 重原子 ≤ 6 — 极简分子，直接 terminal
2. CS score ≤ threshold — 低于配置阈值
3. 无可断键 + 无 FGI — 模板无法处理

目标分子（is_target=true）即使满足 terminal 条件也不会自动跳过，
会进入 standard 流程让 LLM 决策。

---

## 6. 输出文件结构

`export` 命令将结果导出到 `output/YYYYMMDD_HHMMSS_分子名/`:

```
output/20260220_235417_Losartan/
├── SYNTHESIS_REPORT.html   # 自包含 HTML 可视化报告（核心输出）
├── SYNTHESIS_REPORT.md     # Markdown 报告（带图像引用）
├── report.txt              # 正向合成报告（纯文本）
├── tree.json               # 完整合成树 JSON
├── tree.txt                # 合成树文本渲染
├── terminals.json          # 起始原料清单
├── visualization.json      # nodes/edges 图数据（供前端）
├── session.json            # 完整会话快照（可恢复）
└── images/                 # 分子/反应/合成树 PNG 图像
    ├── mol_0.png
    ├── rxn_1_reaction.png
    └── synthesis_tree.png
```

HTML 报告特点：所有分子结构图和反应图内嵌为 base64，无需外部依赖，单文件可分享。
