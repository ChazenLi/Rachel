# Rachel — 多步逆合成规划引擎

LLM 驱动的计算机辅助逆合成分析系统。通过 JSON 接口让 LLM 与化学工具层交互，
逐步拆解目标分子为简单可得的起始原料。

## 架构

```
Rachel/
├── main/                        # 编排引擎（LLM 交互层）
│   ├── retro_cmd.py             # ★ LLM 命令接口（23 个命令）
│   ├── retro_session.py         # JSON 会话持久化 + 沙盒管理
│   ├── retro_orchestrator.py    # BFS 编排器 + 沙盒 + 终止判定
│   ├── retro_tree.py            # 合成树数据模型（双类型节点图）
│   ├── retro_state.py           # 审计状态（决策记录 + 保护基追踪）
│   ├── retro_report.py          # 正向报告 + 可视化数据
│   ├── retro_output.py          # 结果导出（output/ 目录）
│   └── retro_visualizer.py      # HTML/MD 可视化报告 + 分子/反应图像
│
├── chem_tools/                  # 化学工具层（26 个函数，8 个模块）
│   ├── _rdkit_utils.py          # M0: RDKit 公共工具（内部）
│   ├── mol_info.py              # M1: 分子分析
│   ├── fg_detect.py             # M2: 官能团识别
│   ├── template_scan.py         # M3: 模板扫描
│   ├── bond_break.py            # M4: 断键执行
│   ├── forward_validate.py      # M5: 正向验证
│   ├── cs_score.py              # M6: 复杂度评分
│   ├── fg_warnings.py           # M7: 官能团冲突警告
│   ├── smart_cap.py             # M8: 智能断键推理
│   └── templates/               # 反应模板库（12 个 JSON）
│
├── tools/                       # 独立工具脚本
│   ├── llm_bond_test.py         # LLM 交互断键测试平台
│   ├── llm_retro_platform.py    # LLM 逆合成评估平台
│   ├── retro_accuracy.py        # 逆合成准确率评估
│   ├── run_validation.py        # 批量验证运行器
│   ├── multistep_retro.py       # 多步逆合成工具
│   ├── visualize_reaction.py    # 反应可视化 HTML 生成
│   └── _near_miss.py            # Near-miss 分析工具
│
├── tests/                       # 测试套件（210+ 测试）
│   ├── test_m0-m7_unit.py       # 各模块单元测试
│   ├── test_m0-m7_properties.py # 各模块属性测试
│   ├── test_cross_module_properties.py  # 跨模块属性测试
│   ├── test_e2e_integration.py  # 端到端集成测试
│   └── data_driven/             # 数据驱动验证框架
│
├── data/                        # 数据集
│   ├── USPTO50K/                # USPTO50K 反应数据（11 个分片）
│   └── NEW/                     # 新反应数据
│
├── skill.md                     # LLM 操作手册（快速上手）
├── readme.md                    # 本文件
└── refs.md                      # 技术参考（数据结构 + API 详情）
```

## 快速开始

```python
from Rachel.main import RetroCmd

cmd = RetroCmd("my_session.json")

# 1. 创建会话
cmd.execute("init", {
    "target": "CC(=O)Nc1ccc(O)cc1",
    "name": "Paracetamol",
    "terminal_cs_threshold": 1.5,
})

# 2. 获取第一个分子的上下文
ctx = cmd.execute("next")

# 3. 沙盒试探
cmd.execute("try_precursors", {
    "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
    "reaction_type": "Schotten-Baumann acylation",
})

# 4. 提交决策
cmd.execute("select", {"idx": 0})
cmd.execute("commit", {
    "idx": 0,
    "reasoning": "酰氯法，条件温和，前体简单",
    "confidence": "high",
})

# 5. 继续直到完成
ctx = cmd.execute("next")  # → queue_empty
cmd.execute("finalize", {"summary": "一步合成完成"})
cmd.execute("export", {"name": "Paracetamol"})
```

## 核心设计

- 双类型节点图 — MoleculeNode + ReactionNode，canonical SMILES 去重
- BFS 编排 — 广度优先展开，自动处理 terminal 分子
- 沙盒机制 — try_* 不写入树，commit 后才持久化
- 分层上下文 — compact/full/status/tree 四级，控制 token 消耗
- JSON 持久化 — 单文件保存完整状态，对话断了可恢复
- LLM 自提前体 — 不依赖模板，LLM 可基于化学知识直接提出前体
- 智能断键推理 — smart_cap 基于化学环境自动推断 capping 方案
- 可视化报告 — 自包含 HTML 报告，分子/反应图像内嵌 base64

## 依赖

- Python 3.10+
- RDKit
- numpy
- Pillow（可视化报告图像生成）

## 测试

```bash
# 运行全部单元/属性测试（210+ 测试）
pytest Rachel/tests/ -m "not data_driven"

# 运行数据驱动验证
pytest Rachel/tests/data_driven/ --dataset=Rachel/data/USPTO50K/datasetBTF1.csv --sample-size=100 -m data_driven
```

## 文档

- `skill.md` — LLM 操作手册，命令速查，工作流指南
- `refs.md` — 技术参考，数据结构定义，API 详情
- `chem_tools/README.md` — 化学工具层详细设计文档
