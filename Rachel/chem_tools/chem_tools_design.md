# Rachel 化学工具层重构设计方案

## 一、设计目标

将原 `chem_analyze` 的 30+ 文件、数百个函数精简合并为 **7 个核心模块**，消除重复代码，
统一 RDKit 调用路径，保留全部模板数据（templates/），为 LLM 逆合成决策提供纯工具层支撑。

**核心原则**：
- 工具层只做计算和事实输出，不做决策提示
- 能用 RDKit 原生 API 的不再封装
- 每个模块职责单一，模块间通过数据字典通信
- 所有评分标注 `is_heuristic: true`，不伪装精确度

---

## 二、模块架构总览

```
chem_tools/
├── __init__.py                    # 公开 API 导出
├── mol_info.py                    # M1: 分子信息与图结构
├── fg_detect.py                   # M2: 官能团识别与匹配
├── template_scan.py               # M3: 反应模板扫描与匹配
├── bond_break.py                  # M4: 断键执行与前体生成
├── forward_validate.py            # M5: 正向验证与原子守恒
├── cs_score.py                    # M6: CS 合成难度评分
├── fg_warnings.py                 # M7: 官能团排斥与保护警告
├── _rdkit_utils.py                # 内部: RDKit 公共工具函数
└── templates/                     # 数据层: 原样保留全部 JSON 模板
    ├── reactions.json
    ├── functional_groups.json
    ├── fg_compatibility_matrix.json
    ├── protecting_groups.json
    ├── selectivity_conflicts.json
    ├── byproducts.json
    ├── reaction_role_requirements.json
    ├── reactive_sites.json
    ├── structural_alerts.json
    ├── dangerous_combos.json
    ├── known_scaffolds.json
    └── selectivity_reactivity.json
```


**模块依赖关系**（单向，无循环）：

```
_rdkit_utils.py  ← 所有模块共用
      ↑
mol_info.py      ← fg_detect, template_scan, cs_score
      ↑
fg_detect.py     ← template_scan, fg_warnings, forward_validate
      ↑
template_scan.py ← bond_break, forward_validate
      ↑
bond_break.py    ← (被 LLM 调用)
forward_validate.py ← (被 LLM 调用)
cs_score.py      ← (被 LLM 调用)
fg_warnings.py   ← (被 LLM 调用)
```

**原模块 → 新模块映射**：

| 原模块 | 合并到 |
|--------|--------|
| molecule_property_analyzer.py | mol_info + fg_detect + cs_score |
| rdkit_analysis.py | mol_info |
| molecular_graph_engine.py | mol_info |
| graph_constants.py | _rdkit_utils |
| rdkit_operations.py | bond_break |
| disconnection_proposer.py | template_scan |
| disconnection_templates.py | template_scan |
| reaction_pattern_analyzer.py | template_scan + fg_detect |
| retro_template_executor.py | bond_break |
| graph_based_proposer.py | 删除（功能分散到 template_scan + bond_break） |
| precursor_builder.py | bond_break |
| fused_ring_retro.py | bond_break |
| forward_synthesis_validator.py | forward_validate |
| forward_template_executor.py | forward_validate |
| validation_engine.py | forward_validate |
| validation_dimensions.py | forward_validate |
| validation_bond_topology.py | forward_validate |
| atom_mapping_analyzer.py | forward_validate（可选依赖 rxnmapper） |
| synthesis_accessibility_scorer.py | cs_score |
| scaffold_analyzer.py | cs_score |
| selectivity_risk_analyzer.py | fg_warnings |
| protecting_group_detector.py | fg_warnings |
| known_scaffold_matcher.py | mol_info |
| ring_disconnection_analyzer.py | template_scan |
| visual_diagnostic.py | 删除（图像渲染不属于工具层核心） |

---

## 三、各模块详细设计


### M0: `_rdkit_utils.py` — RDKit 公共工具

**职责**：所有模块共用的 RDKit 基础操作，消除各模块重复的 `_parse`、`_canonical`、`_tanimoto` 等。

```python
# === 公开函数 ===

def parse_mol(smiles: str) -> Optional[Chem.Mol]
    """SMILES → RDKit Mol，带 Sanitize，失败返回 None。内置 LRU 缓存。"""

def canonical(smiles: str) -> Optional[str]
    """SMILES → canonical SMILES。"""

def validate_smiles(smiles: str) -> Tuple[bool, str]
    """验证 SMILES 合法性（长度、字符、嵌套深度、RDKit 解析）。"""

def smarts_match(mol: Chem.Mol, smarts: str) -> List[Tuple[int, ...]]
    """SMARTS 子结构匹配，返回所有匹配元组。内置 SMARTS 编译缓存。"""

def tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float
    """Morgan fingerprint (r=2, 2048-bit) Tanimoto 相似度。"""

def mol_formula_counter(mol: Chem.Mol) -> Counter
    """原子计数（含隐式 H），返回 Counter({'C': 6, 'H': 12, ...})。"""

def load_template(filename: str) -> dict
    """加载 templates/ 下的 JSON 文件，过滤 __comment 键。"""

# === 常量 ===
ELECTRONEGATIVITY: Dict[str, float]     # 元素电负性
VALENCE_ELECTRONS: Dict[str, int]       # 价电子数
FRAGMENT_TEMPLATES: Dict[str, str]      # 常用片段 SMILES（"OH", "NH2", "Br" 等）
```

**合并来源**：graph_constants.py 的常量 + 各模块重复的 `_parse`/`_canonical`/`_tanimoto`/`load_json`

---


### M1: `mol_info.py` — 分子信息与图结构

**职责**：分析分子的静态属性，提供原子/键/环/立体化学信息，生成分子图表示。
LLM 用这些信息理解分子结构、选择断键位点。

```python
# === 公开函数 ===

def analyze_molecule(smiles: str) -> Dict[str, Any]
    """
    主入口：全面分析分子属性。

    输入: str SMILES
    输出: {
        "smiles": str,              # canonical SMILES
        "formula": str,             # 分子式
        "mw": float,               # 分子量
        "descriptors": {            # RDKit 描述符
            "logP", "HBD", "HBA", "TPSA",
            "rotatable_bonds", "Fsp3",
            "ring_count", "aromatic_ring_count", "heterocycle_count"
        },
        "atoms": [{                 # 每个原子的信息
            "idx": int, "element": str, "charge": int,
            "hybridization": str, "in_ring": bool,
            "aromatic": bool, "num_hs": int,
            "neighbors": [int, ...]
        }, ...],
        "bonds": [{                 # 每个键的信息
            "idx": int, "atoms": [int, int],
            "bond_type": str,       # SINGLE/DOUBLE/TRIPLE/AROMATIC
            "in_ring": bool, "aromatic": bool,
            "conjugated": bool
        }, ...],
        "rings": [{                 # 环信息
            "atoms": [int, ...], "size": int,
            "aromatic": bool, "heterocyclic": bool
        }, ...],
        "stereo": {                 # 立体化学
            "chiral_centers": [{"idx": int, "cfg": str}, ...],
            "double_bond_stereo": [{"atoms": [int,int], "stereo": str}, ...]
        },
        "symmetry": {               # 对称性
            "ratio": float,         # unique_ranks / n_heavy, 越低越对称
            "level": str            # none/low/moderate/high
        },
        "scaffold": {               # 骨架拓扑
            "murcko": str,          # Murcko 骨架 SMILES
            "tags": [str, ...],     # fused_bicyclic/bridged/spiro/macrocycle 等
            "fused_pairs": int,
            "bridged_pairs": int,
            "spiro_pairs": int
        }
    }

    实现要点:
    - 原子/键信息直接遍历 RDKit Mol 对象获取，不经过 NetworkX
    - 环信息用 mol.GetRingInfo()
    - 立体化学用 Chem.FindMolChiralCenters + FindPotentialStereo
    - 对称性用 Chem.CanonicalRankAtoms(breakTies=False)
    - 骨架用 MurckoScaffold.GetScaffoldForMol
    """

def get_atom_info(mol: Chem.Mol, idx: int) -> Dict[str, Any]
    """单个原子的详细信息。直接调用 RDKit Atom API。"""

def get_bond_info(mol: Chem.Mol, i: int, j: int) -> Dict[str, Any]
    """两个原子间键的详细信息。直接调用 mol.GetBondBetweenAtoms。"""

def match_known_scaffolds(smiles: str, max_hits: int = 5) -> Dict[str, Any]
    """
    与已知骨架库匹配。

    输入: str SMILES
    输出: {
        "ok": bool,
        "hits": [{
            "scaffold_id": str, "name": str, "class": str,
            "method": str,      # smarts/murcko_exact/murcko_substructure
            "coverage": float,  # 匹配原子覆盖率
            "confidence": float,
            "matched_atoms": [int, ...],
            "reference_route_tags": [str, ...]
        }, ...]
    }

    实现要点:
    - 加载 templates/known_scaffolds.json
    - 先 SMARTS 匹配，再 Murcko 回退
    - confidence 标注 is_heuristic
    """
```

**合并来源**：rdkit_analysis.py + molecular_graph_engine.py（查询部分）+ scaffold_analyzer.py + known_scaffold_matcher.py + molecule_property_analyzer.py（基础属性部分）

**删除的冗余**：
- NetworkX 图层（已在旧代码中废弃）
- `get_adjacency_matrix`（LLM 不需要邻接矩阵，bonds 列表等价）
- `get_llm_prompt` / `get_llm_context`（工具层不生成 prompt）
- `_compact_output` / `_build_synthesis_hint`（LLM 自行提取信息）
- druglikeness 判定（不属于逆合成核心）

---


### M2: `fg_detect.py` — 官能团识别与匹配

**职责**：识别分子中的所有官能团，提供官能团 → 可能反应的映射。
这是 template_scan 和 fg_warnings 的基础数据源。

```python
# === 公开函数 ===

def detect_functional_groups(smiles: str) -> Dict[str, Any]
    """
    识别分子中的所有官能团。

    输入: str SMILES
    输出: {
        "ok": bool,
        "smiles": str,
        "groups": {
            "aldehyde": {"count": 1, "atoms": [[3]]},
            "primary_alcohol": {"count": 2, "atoms": [[7], [12]]},
            ...
        }
    }

    实现要点:
    - 加载 templates/functional_groups.json
    - 对每个 SMARTS 执行 mol.GetSubstructMatches
    - 只输出命中的官能团（零命中不输出）
    - bridgehead_atom 按中心原子去重
    """

def detect_reactive_sites(smiles: str) -> Dict[str, Any]
    """
    识别反应活性位点（亲电/亲核/还原/偶联等）。

    输入: str SMILES
    输出: {
        "ok": bool,
        "sites": {
            "electrophilic": {
                "electrophilic_C_ketone": {"count": 1, "atoms": [[3,4]]},
                ...
            },
            "cross_coupling": { ... },
            ...
        }
    }

    实现要点:
    - 加载 templates/reactive_sites.json
    - 按类别前缀分组（electrophilic_/basic_/radical_/...）
    """

def detect_protecting_groups(smiles: str) -> Dict[str, Any]
    """
    检测已存在的保护基。

    输入: str SMILES
    输出: {
        "status": "ok" | "no_pg_found" | "invalid_smiles",
        "detected": [{
            "name": "TBS", "type": "hydroxyl",
            "count": 1, "removal": "TBAF/THF",
            "orthogonal_to": ["Boc", "Fmoc"]
        }, ...],
        "summary": {
            "has_protected_amine": bool,
            "has_protected_hydroxyl": bool,
            ...
            "total_pg_count": int
        },
        "orthogonality_notes": [str, ...]
    }

    实现要点:
    - 加载 templates/protecting_groups.json
    - 预编译 SMARTS，模块加载时执行一次
    - 子结构包含去重（TBS 命中时移除 TMS）
    - 双向正交性检查
    """

def get_fg_reaction_mapping(groups: Dict) -> Dict[str, List[str]]
    """
    官能团 → 可能反应类型的映射。

    输入: detect_functional_groups 的 groups 输出
    输出: {
        "aldehyde": ["nucleophilic_addition", "wittig_olefination",
                     "reductive_amination", "aldol_condensation"],
        "primary_alcohol": ["oxidation_to_aldehyde", "esterification",
                           "mesylation", "tosylation"],
        ...
    }

    实现要点:
    - 内置映射表（从原 _TRANSFORMATION_LABELS 精简而来）
    - 只返回分子中实际存在的官能团对应的反应
    """
```

**合并来源**：molecule_property_analyzer.py（官能团检测部分）+ protecting_group_detector.py + reaction_pattern_analyzer.py（活性位点部分）+ selectivity_risk_analyzer.py（官能团检测部分）

**删除的冗余**：
- `_REACTIVE_CATEGORY_MAP`（分类直接用前缀）
- `_ELECTROPHILIC_FG_SET`（合并到 get_fg_reaction_mapping）
- 危险组合检测移到 fg_warnings

---


### M3: `template_scan.py` — 反应模板扫描与匹配

**职责**：扫描 reactions.json 中的 680 条模板，找出目标分子可能适用的反应，
识别可断键位点及对应的逆合成模板。LLM 基于此输出做断键决策。

```python
# === 数据类 ===

@dataclass
class TemplateMatch:
    template_id: str            # 如 "C08_suzuki_retro"
    name: str                   # 如 "Suzuki Coupling (Retro)"
    category: str               # 如 "cc_formation"
    type: str                   # "retro" | "forward"
    prereq_smarts: str          # 前提条件 SMARTS
    rxn_smarts: str             # 反应 SMARTS
    matched_atoms: List[Tuple[int, ...]]  # prereq 在目标分子上的匹配
    confidence: float           # 模板置信度 (is_heuristic)
    broken_bonds: List[Tuple[int, int]]   # 该模板实际断开的键（映射到目标分子）
    incompatible_groups: List[str]        # 不兼容官能团 SMARTS
    conditions: Optional[str]   # 反应条件描述

# === 公开函数 ===

def scan_applicable_reactions(
    smiles: str,
    mode: str = "retro",       # "retro" | "forward" | "both"
) -> Dict[str, Any]
    """
    扫描目标分子适用的反应模板。

    输入: str SMILES, str mode
    输出: {
        "ok": bool,
        "smiles": str,
        "matches": [TemplateMatch, ...],  # 按 confidence 降序
        "by_category": {                  # 按反应类别分组
            "cc_formation": [TemplateMatch, ...],
            "oxidation_reduction": [...],
            ...
        },
        "summary": {
            "total_matches": int,
            "categories": [str, ...]
        }
    }

    实现要点:
    - 遍历 reactions.json，按 mode 过滤 type
    - 对每个模板：prereq SMARTS 匹配 → 命中则加入结果
    - 从 rxn SMARTS 提取 broken_bonds（反应中心键）
    - confidence 基于匹配质量（prereq 原子覆盖率）
    """

def find_disconnectable_bonds(smiles: str) -> Dict[str, Any]
    """
    识别可断键位点，每个键关联匹配的逆合成模板。

    输入: str SMILES
    输出: {
        "ok": bool,
        "smiles": str,
        "bonds": [{
            "atoms": [int, int],
            "bond_type": str,
            "in_ring": bool,
            "templates": [{             # 该键匹配的逆合成模板
                "template_id": str,
                "name": str,
                "category": str,
                "confidence": float
            }, ...],
            "heuristic_score": float    # 综合断键推荐分
        }, ...],
        "unmatched_bonds": [{           # 有战略价值但无模板匹配的键
            "atoms": [int, int],
            "reason": str
        }, ...]
    }

    实现要点:
    - 调用 scan_applicable_reactions(mode="retro")
    - 对每个模板，从 rxn SMARTS 提取 broken map pairs
    - 通过 prereq 子结构匹配映射到目标分子的真实键
    - 按键合并：一个键下挂多个匹配模板
    - heuristic_score = sqrt(bond_priority * best_template_confidence)
    - 战略键识别：FGI 邻近键、杂原子-碳键、环外键等（Corey 启发式）
    """

def extract_broken_bonds(rxn_smarts: str, prereq_smarts: str,
                         mol: Chem.Mol) -> Set[Tuple[int, int]]
    """
    从 rxn SMARTS 提取反应中心键，映射到目标分子的真实原子索引。

    输入: rxn SMARTS, prereq SMARTS, 目标分子 Mol
    输出: set of (atom_i, atom_j) 元组

    实现要点:
    - 解析 rxn SMARTS 的 reactant/product 模板
    - 比较两侧 atom map number 对之间的键拓扑
    - reactant 有键但 product 无键 → 断键
    - 通过 prereq 子结构匹配映射到真实原子索引
    """
```

**合并来源**：disconnection_proposer.py + disconnection_templates.py + reaction_pattern_analyzer.py（模板匹配部分）+ retro_template_executor.py（模板解析部分）+ ring_disconnection_analyzer.py + graph_based_proposer.py（上下文生成部分）

**删除的冗余**：
- `GraphBasedProposer` 类（拆分为 find_disconnectable_bonds + bond_break 的执行函数）
- `GraphProposalContext` / `ProposalResult` 数据类（用普通 Dict 替代）
- `TemplateMatcher` 类（合并到 scan_applicable_reactions）
- `get_llm_prompt` / `get_json_context`（工具层不生成 prompt）
- `_bond_overlaps_prereq`（已废弃的旧匹配逻辑）
- `ReactionCategory` 枚举（直接用字符串）

---


### M4: `bond_break.py` — 断键执行与前体生成

**职责**：接收 LLM 的断键决策，在 RDKit Mol 上执行操作，生成前体 SMILES。
这是整个系统中唯一修改分子结构的模块。

```python
# === 数据类 ===

@dataclass
class BreakResult:
    success: bool
    precursors: List[str]       # 前体 SMILES 列表
    operation_log: List[str]    # 操作日志
    error: Optional[str] = None

# === 公开函数 ===

def execute_disconnection(
    smiles: str,
    bond: Tuple[int, int],
    reaction_type: str,
    custom_fragments: Optional[List[Tuple[str, str]]] = None,
) -> BreakResult
    """
    执行单键断裂，生成前体。

    输入:
        smiles: 目标分子 SMILES
        bond: (atom_i, atom_j) 要断开的键
        reaction_type: 反应类型（template_id 或自由文本）
        custom_fragments: 可选，自定义片段 [(fragment_smiles, target_side), ...]
            target_side: "i" 接到 atom_i, "j" 接到 atom_j

    输出: BreakResult

    执行优先级:
    1. custom_fragments（LLM 显式指定片段）
    2. rxn SMARTS 模板执行（从 reactions.json 查找匹配模板，
       用 AllChem.ReactionFromSmarts 执行逆合成）
    3. 简单断键 + 加 H（兜底：断开键，两端各补一个 H）

    实现要点:
    - 在 RWMol 上操作：RemoveBond → GetMolFrags → 各片段 SanitizeMol
    - rxn SMARTS 执行：AllChem.ReactionFromSmarts → RunReactants
    - 片段添加：CombineMols + AddBond
    """

def execute_operations(
    smiles: str,
    operations: List[Dict[str, Any]],
    reaction_category: Optional[str] = None,
) -> BreakResult
    """
    执行操作序列（LLM 自主决策模式，支持复杂断键）。

    输入:
        smiles: 目标分子 SMILES
        operations: 操作列表，每个操作:
            {"op": "disconnect", "bond": [i, j]}
            {"op": "add_fragment", "atom": i, "fragment": "Br", "bond_type": "SINGLE"}
            {"op": "change_bond_order", "bond": [i, j], "new_order": 2}
            {"op": "ring_open", "bond": [i, j]}
            {"op": "modify_atom", "atom": i, "charge": -1}
            {"op": "fused_ring_retro", "bonds": [[i,j],[k,l]],
             "caps": {...}, "bond_order_changes": [...]}
        reaction_category: 反应类别（用于角色检查）

    输出: BreakResult

    实现要点:
    - 按顺序在 RWMol 上执行每个操作
    - fused_ring_retro 使用专用路径处理稠环系统
    - 最终 GetMolFrags 拆分片段，SanitizeMol 清理
    - 可选：角色检查（前体是否满足亲核/亲电要求）
    """

def try_retro_template(
    mol: Chem.Mol,
    bond: Tuple[int, int],
    template_id: str,
) -> Optional[BreakResult]
    """
    尝试用指定逆合成模板在指定键上执行。

    输入: Mol, 键, 模板 ID
    输出: BreakResult 或 None（模板不匹配时）

    实现要点:
    - 从 reactions.json 加载模板的 rxn SMARTS
    - AllChem.ReactionFromSmarts 执行
    - 过滤产物：只保留断开指定键的结果
    - SanitizeMol 清理产物
    """

def resolve_template_ids(reaction_name: str) -> List[str]
    """
    根据反应名称模糊匹配模板 ID 列表。

    输入: str 反应名（如 "suzuki", "aldol condensation"）
    输出: [str, ...] 匹配的模板 ID 列表

    实现要点:
    - 名称标准化（小写、去空格、去连字符）
    - 在 reactions.json 的 name 字段中模糊搜索
    """
```

**合并来源**：rdkit_operations.py + retro_template_executor.py + precursor_builder.py + fused_ring_retro.py + graph_based_proposer.py（执行部分）+ molecular_graph_engine.py（execute_disconnection/execute_operations）

**删除的冗余**：
- `MolecularGraphEngine` 类（Facade 层不再需要，直接调用函数）
- `preview_bond_atoms`（LLM 从 mol_info 获取键信息）
- `build_coupling_pair` / `build_ring_open` / `build_from_proposal` / `build_multi_break` / `suggest_and_build`（这些是 precursor_builder 的便捷封装，统一用 execute_disconnection + execute_operations 替代）
- `validate_mol_operations`（合并到 execute_operations 内部）

---


### M5: `forward_validate.py` — 正向验证与原子守恒

**职责**：LLM 提出断键方案和前体后，验证正向反应可行性。
包含原子守恒检查、模板正向执行、骨架对齐、键拓扑验证。

```python
# === 公开函数 ===

def validate_forward(
    precursors: List[str],
    target: str,
    template_id: Optional[str] = None,
    reaction_smarts: Optional[str] = None,
    reaction_category: Optional[str] = None,
    byproduct_smiles: Optional[List[str]] = None,
) -> Dict[str, Any]
    """
    正向合成可行性验证（主入口）。

    输入:
        precursors: 前体 SMILES 列表
        target: 目标分子 SMILES
        template_id: 可选，reactions.json 模板 ID
        reaction_smarts: 可选，直接传入的反应 SMARTS
        reaction_category: 可选，反应类别
        byproduct_smiles: 可选，LLM 声明的副产物 SMILES

    输出: {
        "ok": bool,
        "target": str,
        "precursors": [str, ...],
        "checks": {
            "atom_balance": { ... },        # 原子守恒
            "template_execution": { ... },  # 模板正向执行
            "scaffold_alignment": { ... },  # MCS 骨架对齐
            "bond_topology": { ... },       # 键变化拓扑
            "fg_compatibility": { ... },    # 官能团兼容性
        },
        "assessment": {
            "feasibility_score": float,     # 综合可行性分 [0, 1]
            "pass": bool,                   # 是否通过 hard gate
            "hard_fail_reasons": [str, ...] | None,
            "score_breakdown": {
                "atom_balance": float,
                "template_match": float,
                "scaffold_alignment": float,
                "bond_topology": float,
                "fg_compatibility": float
            }
        }
    }

    验证步骤:
    1. 原子守恒（check_atom_balance）
    2. 模板正向执行（execute_forward_template）
    3. MCS 骨架对齐（check_scaffold_alignment）
    4. 键变化拓扑（check_bond_topology）
    5. 官能团兼容性（check_fg_compatibility）
    6. 聚合评分 + hard gate 判定
    """

def check_atom_balance(
    precursors: List[str],
    target: str,
    byproduct_smiles: Optional[List[str]] = None,
    reaction_category: Optional[str] = None,
) -> Dict[str, Any]
    """
    原子守恒检查（两层小分子损失设计）。

    输入:
        precursors: 前体 SMILES 列表
        target: 目标分子 SMILES
        byproduct_smiles: LLM 声明的副产物（如 ["O", "[Br-]"]）
        reaction_category: 反应类别（用于 Tier 1 损失查表）

    输出: {
        "balanced": bool,           # 两层损失扣除后是否完全守恒
        "balance_score": float,     # 连续评分 [0, 1]，= 1 - unexplained/total_heavy
        "precursor_atoms": Counter, # 前体侧原子计数
        "product_atoms": Counter,   # 产物侧原子计数（target + byproducts）
        "excess": Counter,          # 产物多出的原子（raw）
        "deficit": Counter,         # 前体多出的原子（raw）
        "adjusted_excess": dict,    # 两层损失扣除后仍多出的原子
        "adjusted_deficit": dict,   # 两层损失扣除后仍缺少的原子
        "severe_imbalance": bool,   # 产物侧非 H 多出 > 4（hard gate）
        "skeleton_imbalance": bool, # 产物比前体多出 C/N/S（hard gate）
        "note": str
    }

    两层损失设计:
    - Tier 1: _SMALL_MOL_LOSSES[reaction_category] → 类别特定损失（19 类反应）
    - Tier 2: _UNIVERSAL_LOSS_KEYS → 通用非骨架损失（H₂O, HCl, HBr, HI, HF, H₂）
    - 两层依次应用于 deficit 和 excess 两侧
    - _LOSS_COUNTERS 共 13 种小分子（含 N₂, CO₂, MeOH, EtOH, AcOH, NH₃, SO₂）

    单向不平衡检测:
    - skeleton_imbalance: 仅产物侧多出 C/N/S 时触发
      （前体侧多出是正常的试剂离去基团，如 Wittig 的 Ph₃P=O）
    - severe_imbalance: 仅产物侧非 H 多出 > 4 时触发
      （前体侧多出通过 balance_score 连续评分软惩罚）

    连续 balance_score:
    - 公式: 1.0 - (unexplained_atoms / total_heavy_atoms)
    - 取代旧的二值 balanced/not 判定
    - 直接参与 validate_forward 的加权平均（权重 0.25）
    - 典型值: 完全守恒=1.0, Wittig≈0.77, Suzuki≈0.81, Stille≈0.5
    """

def execute_forward_template(
    precursors: List[str],
    target: str,
    template_id: Optional[str] = None,
    reaction_smarts: Optional[str] = None,
) -> Dict[str, Any]
    """
    用前体跑正向反应 SMARTS，检查目标是否在产物中。

    输入: 前体列表, 目标, 可选模板
    输出: {
        "attempted": bool,
        "target_in_products": bool,
        "tanimoto_to_target": float,
        "products_generated": [str, ...],
        "best_match": str
    }

    实现要点:
    - 从 reactions.json 查找正向模板（retro 模板自动翻转）
    - AllChem.ReactionFromSmarts → RunReactants
    - 产物清理：三级回退 SanitizeMol
    - 与目标比较：精确匹配 + Tanimoto 相似度
    """

def check_scaffold_alignment(
    precursors: List[str],
    target: str,
    reaction_category: Optional[str] = None,
) -> Dict[str, Any]
    """
    MCS 骨架对齐验证。

    输出: {
        "aligned": bool,
        "coverage_ratio": float,    # MCS 覆盖目标原子的比例
        "mcs_level": int,           # 1=普通, 2=缩合, 3=成环
        "warning": Optional[str]
    }

    实现要点:
    - rdFMCS.FindMCS 计算最大公共子结构
    - 根据 reaction_category 调整阈值（成环反应更宽松）
    """

def check_bond_topology(
    precursors: List[str],
    target: str,
    reaction_category: Optional[str] = None,
) -> Dict[str, Any]
    """
    键变化拓扑验证（新键类型、成环可行性、稠环融合）。

    输出: {
        "pass": bool,
        "has_hard_fail": bool,
        "violations": [{"reason": str, "severity": str}, ...],
        "summary": str
    }
    """
```

**合并来源**：validation_engine.py + validation_dimensions.py + validation_bond_topology.py + forward_template_executor.py + atom_mapping_analyzer.py

**关键简化**：
- 原 7 步验证精简为 5 步（移除独立的 SA 改善评估和原子映射验证）
- SA 改善评估移到 cs_score 模块，按需调用
- 原子映射验证（rxnmapper）作为可选增强，不在主流程中
- `_build_llm_context` 删除（工具层不生成 LLM 上下文）
- `_compact_facts` 删除（始终返回完整 checks）
- 权重和阈值集中定义，不散落在各函数中

---


### M6: `cs_score.py` — CS (Complicated Score) 合成难度评分

**职责**：自定义 CS 评分系统，评估分子的合成难度。
用于判断前体是否足够简单可作为合成终点，以及评估合成路线进度。

**CS 评分设计**（替代原 SA Score，更透明可控）：

```python
# === CS 评分维度与权重 ===

CS_DIMENSIONS = {
    "size":       0.55,  # 分子大小（重原子数，非 MW）
    "ring":       0.65,  # 环系统拓扑（稠合/桥环/螺环/大环）
    "stereo":     0.55,  # 立体化学负担（手性中心数、密度、连续性）
    "hetero":     0.40,  # 杂原子密度与多样性
    "symmetry":  -0.20,  # 对称性折扣（高对称 → 减分）
    "fg_density": 0.35,  # 官能团密度（保护基 + 复杂度官能团）
}

# CS Score 范围: [1.0, 10.0]，越低越容易合成
# 理论最大值约 8.75
# 分类阈值:
CS_TRIVIAL = 2.5    # <= 2.5: 合成极简单，可作为终点
CS_MODERATE = 6.0   # <= 6.0: 中等难度
                     # > 6.0: 复杂

# === 公开函数 ===

def compute_cs_score(smiles: str) -> Dict[str, Any]
    """
    计算分子的 CS (Complicated Score) 合成难度评分。

    输入: str SMILES
    输出: {
        "cs_score": float,          # [1.0, 10.0]
        "classification": str,      # "trivial" | "moderate" | "complex"
        "is_terminal": bool,        # cs_score <= CS_TRIVIAL
        "is_heuristic": true,       # 始终标注
        "breakdown": {
            "size": float,          # 分子大小贡献
            "ring": float,          # 环系统贡献
            "stereo": float,        # 立体化学贡献
            "hetero": float,        # 杂原子贡献
            "symmetry": float,      # 对称性折扣（负值）
            "fg_density": float     # 官能团密度贡献
        }
    }

    各维度计算:

    1. size: min(log2(n_heavy/6) * 1.2, 3.0)
       - 基于重原子数（非 MW），避免 CF₃ 等基团虚增
       - n_heavy <= 6 → 0
       - 10 heavy → 0.35, 20 → 1.05, 30 → 1.50, 50 → 2.10

    2. ring:
       - 独立环: min(n * 0.4, 1.6)
       - 稠合: n_fused * 0.8
       - 桥环: n_bridged * 1.5
       - 螺环: n_spiro * 1.0
       - 大环(>8): n_macro * 1.5
       - 杂环微调: min(n_hetero_ring * 0.2, 1.0)
       - 上限 5.0

    3. stereo:
       - 基础: 0.4 * sqrt(n_chiral)
       - 密度 > 0.15: (density - 0.15) * 5.0
       - 连续手性中心: min(n_consecutive * 0.3, 1.5)
       - 上限 3.5

    4. hetero:
       - 密度 > 0.1: (density - 0.1) * 3.0
       - 种类多样性: min(n_types * 0.15, 0.8)
       - 上限 2.5

    5. symmetry (折扣):
       - ratio < 0.8: (0.8 - ratio) * 2.0
       - ratio = unique_ranks / n_heavy

    6. fg_density:
       - 保护基匹配数 * 0.35 (TMS/Boc/Cbz/Ac/Bn/Ts/Moc/缩醛)
       - 复杂度官能团匹配数 * 0.25 (酰胺/磺酰胺/氨基甲酸酯/脲/酯/酰氯/硼酸)
       - 上限 2.5
    """

def classify_complexity(smiles: str) -> Dict[str, Any]
    """
    快速分类：trivial / moderate / complex。

    输入: str SMILES
    输出: {
        "cs_score": float,
        "classification": str,
        "is_terminal": bool
    }
    """

def score_progress(
    intermediate: str,
    target: str,
) -> Dict[str, Any]
    """
    评估中间体到目标的合成进度。

    输入: 中间体 SMILES, 目标 SMILES
    输出: {
        "ok": bool,
        "progress_score": float,    # [0, 1]，越高越接近目标
        "classification": str,      # "close" (>=0.7) / "moderate" (>=0.4) / "early"
        "raw_similarity": float,    # 原始 Tanimoto
        "adjustments": {
            "deprotection_boost": float,    # 虚拟脱保护修正
            "substruct_bonus": float,       # 子结构加分
            "fg_transform_bonus": float     # 官能团转换修正
        },
        "cs_score": float           # 中间体的 CS Score
    }

    实现要点:
    - Morgan fingerprint Tanimoto 相似度
    - 虚拟脱保护：移除常见保护基后重新计算
    - 子结构加分：中间体是目标子结构时
    - 官能团转换修正：已知前体→目标官能团对
    """

def batch_score_progress(
    intermediates: List[str],
    target: str,
) -> List[Dict[str, Any]]
    """批量评估，预计算 target 侧数据避免重复。"""
```

**合并来源**：synthesis_accessibility_scorer.py（全部）+ scaffold_analyzer.py（复杂度部分）+ molecule_property_analyzer.py（SA 部分）

**关键变更**：
- 重命名 SA Score → CS Score (Complicated Score)，更准确反映含义
- 移除 RDKit sascorer 依赖（原 blended 策略），纯启发式
- 评分维度和权重显式定义，可调可追溯
- 始终标注 `is_heuristic: true`

---


### M7: `fg_warnings.py` — 官能团排斥与保护警告

**职责**：检测官能团间的化学冲突、反应条件不兼容、保护基需求。
在 LLM 选择反应路线时提供安全性警告。

```python
# === 公开函数 ===

def check_fg_conflicts(smiles: str) -> Dict[str, Any]
    """
    检测分子内官能团间的化学冲突。

    输入: str SMILES
    输出: {
        "ok": bool,
        "conflicts": [{
            "group_a": str,         # 如 "aldehyde"
            "group_b": str,         # 如 "primary_amine"
            "reaction_type": str,   # 如 "imine_formation"
            "severity": str,        # "high" | "medium" | "low"
            "note": str             # 如 "醛与伯胺自发缩合形成亚胺"
        }, ...],
        "competing_groups": [{
            "category": str,        # 如 "electrophilicity_competition"
            "groups": [str, ...],   # 竞争的官能团
            "note": str
        }, ...],
        "dangerous_combos": [str, ...]  # 危险组合警告
    }

    实现要点:
    - 先调用 fg_detect.detect_functional_groups 获取官能团
    - 加载 templates/selectivity_conflicts.json 检查冲突对
    - 加载 templates/dangerous_combos.json 检查危险组合
    """

def check_reaction_compatibility(
    smiles: str,
    reaction_category: str,
) -> Dict[str, Any]
    """
    检查分子中的官能团是否与指定反应条件兼容。

    输入: str SMILES, str 反应类别（如 "strong_oxidation"）
    输出: {
        "compatible": bool,
        "warnings": [{
            "fg": str,             # 不兼容的官能团
            "level": str,          # "risky" | "forbidden"
            "reaction_family": str,
            "note": str
        }, ...]
    }

    实现要点:
    - 加载 templates/fg_compatibility_matrix.json
    - 检查分子中每个官能团在该反应条件下的兼容性
    """

def suggest_protection_needs(
    smiles: str,
    planned_reaction: str,
) -> Dict[str, Any]
    """
    根据计划反应，建议需要保护的官能团。

    输入: str SMILES, str 计划反应类别
    输出: {
        "needs_protection": bool,
        "suggestions": [{
            "fg": str,              # 需要保护的官能团
            "reason": str,          # 为什么需要保护
            "recommended_pg": [str, ...],  # 推荐的保护基
            "removal_conditions": [str, ...]
        }, ...],
        "existing_pg": [{           # 已有的保护基
            "name": str,
            "type": str,
            "removal": str
        }, ...],
        "orthogonality_notes": [str, ...]
    }

    实现要点:
    - 调用 fg_detect.detect_functional_groups + detect_protecting_groups
    - 用 fg_compatibility_matrix 找出与 planned_reaction 不兼容的 FG
    - 从 protecting_groups.json 推荐合适的保护基
    - 检查已有保护基与新保护基的正交性
    """

def check_deprotection_safety(
    smiles: str,
    pg_to_remove: str,
) -> Dict[str, Any]
    """
    检查脱保护操作对分子中其他官能团的影响。

    输入: str SMILES, str 要脱除的保护基名称
    输出: {
        "safe": bool,
        "risks": [{
            "affected_fg": str,
            "condition_family": str,
            "level": str            # "risky" | "forbidden"
        }, ...]
    }

    实现要点:
    - 从 protecting_groups.json 获取脱除条件
    - 推断脱除条件涉及的反应家族
    - 用 fg_compatibility_matrix 检查其他 FG 的兼容性
    """
```

**合并来源**：selectivity_risk_analyzer.py + protecting_group_detector.py（正交性/安全性部分）+ molecule_property_analyzer.py（_build_synthesis_hint 中的选择性风险部分）+ validation_dimensions.py（_check_precursor_compatibility）

**删除的冗余**：
- `_build_protection_plan`（合并到 suggest_protection_needs）
- `_differentiation_hint`（合并到 check_fg_conflicts）
- `_PG_REMOVAL_TO_FAMILY` / `_PG_NAME_TO_FG_KEYS`（合并到 check_deprotection_safety 内部）

---


## 四、`__init__.py` 公开 API 导出

```python
"""
chem_tools — Rachel 化学工具层
================================
为 LLM 逆合成决策提供纯计算工具支撑。
"""

# M1: 分子信息
from .mol_info import analyze_molecule, match_known_scaffolds

# M2: 官能团识别
from .fg_detect import (
    detect_functional_groups,
    detect_reactive_sites,
    detect_protecting_groups,
    get_fg_reaction_mapping,
)

# M3: 模板扫描
from .template_scan import (
    scan_applicable_reactions,
    find_disconnectable_bonds,
)

# M4: 断键执行
from .bond_break import (
    execute_disconnection,
    execute_operations,
    try_retro_template,
    resolve_template_ids,
)

# M5: 正向验证
from .forward_validate import (
    validate_forward,
    check_atom_balance,
)

# M6: CS 评分
from .cs_score import (
    compute_cs_score,
    classify_complexity,
    score_progress,
    batch_score_progress,
)

# M7: 官能团警告
from .fg_warnings import (
    check_fg_conflicts,
    check_reaction_compatibility,
    suggest_protection_needs,
    check_deprotection_safety,
)
```

**总计 22 个公开函数**（原系统 60+ 个），覆盖全部核心功能。


---

## 五、LLM 调用工作流

以下是 LLM 逆合成决策中调用工具层的典型流程：

```
Step 1: 分析目标分子
  → analyze_molecule(target_smiles)
  → detect_functional_groups(target_smiles)
  → compute_cs_score(target_smiles)

Step 2: 扫描可能反应 & 可断键
  → find_disconnectable_bonds(target_smiles)
  → check_fg_conflicts(target_smiles)

Step 3: LLM 决策断键（工具层不参与决策）

Step 4: 执行断键
  → execute_disconnection(target, bond, reaction_type)
  或 execute_operations(target, operations, category)

Step 5: 验证
  → validate_forward(precursors, target, template_id, byproducts)
  → check_atom_balance(precursors, target, byproducts)
  → check_reaction_compatibility(target, reaction_category)

Step 6: 评估前体
  → compute_cs_score(precursor)  → is_terminal?
  → score_progress(precursor, target)

Step 7: 如需保护
  → suggest_protection_needs(precursor, planned_reaction)
  → check_deprotection_safety(precursor, pg_name)

Step 8: 递归（对非 terminal 前体重复 Step 1-7）
```

---

## 六、函数数量对比

| 维度 | 原系统 | 新系统（当前） |
|------|--------|---------------|
| Python 文件 | 25 | 9（含 M8 smart_cap） |
| 公开函数 | 60+ | 26 |
| dataclass | 8 | 3 (BreakResult, TemplateMatch, CappingProposal) |
| 内部辅助函数 | 100+ | ~35 |
| 模板文件 | 12 | 12（不变） |

注：原设计为 7 模块 22 函数，实际实现中新增了 M8 smart_cap 模块（2 个函数）和 M4 的 execute_fgi、preview_disconnections（2 个函数），共 26 个公开函数。

---

## 七、关键设计决策说明


### 7.1 为什么删除 GraphBasedProposer / MolecularGraphEngine

原系统有两条并行的逆合成路径：
- `disconnection_proposer` → `reaction_pattern_analyzer` → `precursor_builder`
- `GraphBasedProposer` → `MolecularGraphEngine` → `rdkit_operations`

两条路径功能高度重叠，都是"扫描模板 → 匹配键 → 执行断键 → 生成前体"。
新系统统一为 `template_scan` → `bond_break` 单一路径。

`MolecularGraphEngine` 是一个 Facade 类，所有方法都是对 `rdkit_analysis` 的一对一委托，
没有自己的逻辑。新系统直接暴露函数，不需要 Facade。

### 7.2 为什么用 CS Score 替代 SA Score

原系统的 SA Score 有三个问题：
1. 名称误导：SA = Synthesis Accessibility，容易与"商业可得性"混淆
2. 依赖 RDKit sascorer（基于 PubChem 片段频率），对药物分子系统性偏低
3. blended 策略（75% heuristic + 25% sascorer）增加了不透明性

CS Score (Complicated Score) 纯启发式，维度和权重完全透明，
LLM 可以看到每个维度的贡献，做出更好的判断。

### 7.3 为什么删除 visual_diagnostic

图像渲染（PNG/SVG）不属于化学计算工具层的核心职责。
如果需要可视化，应该在上层（UI 层或 LLM 工具调用层）单独实现。
这样化学工具层可以在无 GUI 环境下运行。

### 7.4 为什么保留全部 12 个模板文件

模板文件是领域知识的结晶（680 条反应 SMARTS、200+ 官能团定义、保护基正交性映射等），
不应该硬编码到 Python 代码中。保留 JSON 格式便于：
- 独立维护和扩展（不需要改 Python 代码）
- 版本控制和 diff
- 其他工具/语言复用

### 7.5 原子映射验证（rxnmapper）的定位

原系统将 atom_mapping_analyzer 作为验证流程的必经步骤，
但 rxnmapper 是可选依赖（需要额外安装），且对某些反应类型置信度偏低。

新系统将其定位为 forward_validate 的可选增强：
- rxnmapper 可用时，作为额外的验证维度
- 不可用时，其他 4 个验证步骤仍然有效
- 不影响 hard gate 判定（除非置信度极低）

### 7.6 断键执行的三级回退策略

`execute_disconnection` 保留原系统的三级优先级，这是经过验证的可靠策略：
1. custom_fragments（LLM 显式指定，最精确）
2. rxn SMARTS 模板执行（reactions.json 680 条模板，覆盖面广）
3. 简单断键 + 加 H（兜底，确保总能返回结果）

### 7.7 forward_validate 的 hard gate 设计

保留原系统的 hard gate 机制（任一 gate 失败 → 整体 fail），
但精简 gate 数量：

| Gate | 条件 | 含义 |
|------|------|------|
| 1 | severe_imbalance | 产物侧非 H 原子多出 > 4 |
| 2 | skeleton_imbalance | 产物比前体多出骨架原子(C/N/S) |
| 3 | scaffold not aligned | MCS 骨架对齐不足 |
| 4 | bond topology violation | 键变化拓扑不可行 |
| 5 | forbidden FG | 含有反应禁忌官能团 |

注：severe_imbalance 和 skeleton_imbalance 均为单向检测（仅产物侧多出时触发）。前体侧多出骨架原子是正常的（试剂离去基团，如 Wittig 的 Ph₃P=O、Stille 的 Bu₃SnCl），通过连续 `balance_score` 软惩罚。

---

## 八、迁移策略

1. 先实现 `_rdkit_utils.py`（公共基础）
2. 按依赖顺序实现：mol_info → fg_detect → template_scan → bond_break → forward_validate → cs_score → fg_warnings
3. 每个模块实现后，用原系统的测试用例验证输出一致性
4. 最后实现 `__init__.py` 导出，替换原 `chem_analyze` 的 import 路径
