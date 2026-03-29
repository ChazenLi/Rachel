<div align="right">

[English](./README.md) | [简体中文](./README.zh-CN.md)

</div>

<div align="center">

<a id="top"></a>

# Rachel

**面向多步逆合成的形式化化学推理框架**

<img alt="Python 3.10+" src="https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white">
<img alt="Active Research" src="https://img.shields.io/badge/Status-活跃研究中-2D6A4F">
<img alt="Multi-Step Retrosynthesis" src="https://img.shields.io/badge/Domain-多步逆合成-8C564B">
<img alt="Workflow" src="https://img.shields.io/badge/Workflow-状态--动作--验证--提交-7B61FF">
<img alt="LLM Strategy Layer" src="https://img.shields.io/badge/LLM-策略层-6F42C1">

<p>
  <a href="#trace-demo-zh">流程追踪演示</a> |
  <a href="#end-to-end-example-zh">端到端示例</a> |
  <a href="#highlights-zh">核心亮点</a> |
  <a href="#selected-molecules-zh">代表性分子</a> |
  <a href="#minimal-quickstart-zh">快速开始</a>
</p>

https://github.com/user-attachments/assets/4dc9990f-00b2-40d8-a8c3-181c6f0c568b

</div>

Rachel 将逆合成规划建模为一个结构化的 `state -> action -> validation -> commit` 过程，而不是一次性的文本生成任务。当前仓库是一个仍在持续使用中的活跃研究代码库，正在为面向 arXiv 的展示做整理与清理。

<a id="trace-demo-zh"></a>
## 流程追踪演示

上方的 trace 是对 Rachel 规划流程的可视化演示，展示了它如何从结构化状态逐步走向已提交的逆合成路线。

<img width="1560" height="1120" alt="trace_final" src="https://github.com/user-attachments/assets/0eca73f1-25c9-4816-b7da-6bbfc24853e3" />

- 这份 trace 的重点不仅是最终路线结果，更是规划行为本身
- 它是理解 Rachel 如何从上下文走向候选探索、验证与提交的最快方式

<a id="end-to-end-example-zh"></a>
## 端到端示例

下图展示了一个完整的路线级对比：PaRoutes 的 ground-truth 方案与 Rachel 在 `n1_366` 这一案例上的生成结果。

<img width="2500" height="4459" alt="n1_366_groundtruth_vs_rachel_annotated_case_en" src="https://github.com/user-attachments/assets/38952d7e-8dc4-4f92-b13c-eee61175b0ec" />

这个示例作为系统层面的定性参考。重点不只是单步反应在化学上是否合理，而是整条路线在完整规划层面上是否保持可解释性与结构一致性。

<a id="highlights-zh"></a>
## 核心亮点

| 能力 | 在 Rachel 中的体现 |
| --- | --- |
| 有状态规划 | Rachel 基于持久化会话状态进行推理，而不是孤立的一次性回答。 |
| 受化学约束的操作空间 | 键断裂（bond disconnection）与 FGI 被视为互补的规划操作。 |
| 提交前沙盒试验 | 候选步骤会先在本地沙盒中尝试，再决定是否写入主路线树。 |
| 受验证门控的执行 | 可行性检查、原子守恒等验证器帮助控制是否提交。 |
| 结构化路线记忆 | 被接受的步骤会成为显式的路线树对象，而不只是自由文本。 |
| LLM 作为策略层 | LLM 负责组织搜索与决策，而不是充当不受约束的化学“黑箱”。 |

## Rachel 是什么

Rachel 结合了以下几个部分：

- 基于会话驱动的规划流程
- 具备化学基础约束的操作器，如键断裂与 FGI
- 提交前的本地沙盒试验
- 由验证器门控的路线构建
- 路线树记忆与可审计状态

目标不仅是提出一条路线，更是让路线构建过程变得可检查、可恢复、可做化学一致性验证。

## 核心工作流

```mermaid
flowchart LR
    A["当前状态"] --> B["候选动作"]
    B --> C["沙盒试验"]
    C --> D["验证"]
    D --> E["提交或拒绝"]
    E --> F["更新后的路线树"]
```

候选动作会先被探索，再写入路线。通过验证的步骤会进入主路线树；被拒绝的尝试也不会简单消失，而是作为有信息价值的规划痕迹被保留下来。

<a id="selected-molecules-zh"></a>
## 代表性分子

Rachel 当前展示了三个定性示例，用于覆盖不同但互补的能力侧面。

<table>
  <tr>
    <td align="center" width="33%">
      <img src="https://github.com/user-attachments/assets/61f7e78b-053c-4ac4-a349-b22c9e5b1ae3" alt="QNTR" width="220"><br>
      <strong>QNTR</strong>
    </td>
    <td align="center" width="33%">
      <img src="https://github.com/user-attachments/assets/e27005c7-9ba1-470b-a038-41d2190e3c72" alt="Losartan" width="220"><br>
      <strong>Losartan</strong>
    </td>
    <td align="center" width="33%">
      <img src="https://github.com/user-attachments/assets/ff2abe54-20c4-427a-8363-b9b6b8634a23" alt="Rivaroxaban" width="220"><br>
      <strong>Rivaroxaban</strong>
    </td>
  </tr>
</table>

| 分子 | 角色 | 路线深度 | 体现的特点 |
| --- | --- | ---: | --- |
| `QNTR` | 具有实验基础的示例 | 2 步 | 一条简短且可解释、并与真实合成经验相关的路线 |
| `Losartan` | 经典药物化学目标 | 4 步 | 体现具有辨识度的药化断裂逻辑与汇聚式路线设计 |
| `Rivaroxaban` | 更深层的类药分子示例 | 5 步 | 展示更长程规划能力与更丰富的转化类型 |

### QNTR

一个有实验基础支撑、并且与我本人合成经验相关的分子。这是一个真正的、针对 **NO2 还原型探针** 的荧光探针，而我已经完成了整个化学合成过程，与条件优化。我们都选择将目标分子拆分为三个主要部分，并采用相同的末端分子、相同的中间体和反应。

尽管 Rachel 的第一个版本存在错误，做出的官能团互变（FGI）处理也不够理想，而且在没有参考文献的现实情况下，环合/开环本身也并不容易预测，更不用说在实验前就跑出合成路线图了。但整个路线和设计策略依然让我感到兴奋。这也确实是一个非常重要的原因，使我开始想把它做成一个在化学上更可行的系统。

#### 我的路线

<img width="954" height="538" alt="1878e5777ea5c79edce765660331f35d" src="https://github.com/user-attachments/assets/1bc9ec91-4137-4fa8-9a42-df25d2af2c0f" />  

#### Rachel 的早期版本路线

<img width="2260" height="2150" alt="synthesis_tree - 副本" src="https://github.com/user-attachments/assets/b03cfe51-d94e-4271-8000-ed0b1712810c" />

#### Rachel 的路线

<img width="2580" height="2150" alt="synthesis_tree" src="https://github.com/user-attachments/assets/00d1a70f-36f9-4077-b0fa-59cea407faac" />


- 2 步路线，起始于 3 个原料
- 适合作为一个紧凑且可解释的规划演示
- 尤其有价值，因为它展示的已接受路线并不只是模板匹配，而体现了策略选择

### Losartan

一个经典的药物化学目标，具有易于识别的汇聚式路线。

- 4 步路线，起始于 4 个原料
- 突出展示 tetrazole formation、N-alkylation 与 Suzuki coupling 等逻辑
- 适合作为许多读者都能快速理解的 benchmark 风格示例

### Rivaroxaban

一个更深层的类药分子示例，具有更丰富的转化组合。

- 5 步路线，起始于 4 个原料
- 突出展示 Buchwald-Hartwig amination、FGI、环化以及酰胺形成
- 有助于说明 Rachel 并不局限于短路线或玩具级案例

### 双药物案例对比

下图将 Losartan 与 Rivaroxaban 两个示例放在同一张带注释的对比图中，为 README 提供更完整的定性展示，帮助读者理解 Rachel 在两类可识别药物样目标上的行为差异。

<img width="3000" height="3755" alt="rivaroxaban_losartan_dual_annotated_en" src="https://github.com/user-attachments/assets/8cb6c479-3f63-41fc-921f-62a565909dd1" />

综合来看，这两个案例形成了清晰对比：

- `Losartan` 强调经典的汇聚式药物化学逻辑
- `Rivaroxaban` 强调更深的路线深度与更丰富的操作器多样性
- 二者组合有助于读者比较路线风格，而不仅是孤立结果

<a id="minimal-quickstart-zh"></a>
## 最小快速开始

当前本地运行默认你已经准备好了主要研究依赖环境，包括 Python 3.10+、RDKit、`numpy` 和 `Pillow`。

```python
from Rachel.main import RetroCmd

cmd = RetroCmd("my_session.json")

cmd.execute(
    "init",
    {
        "target": "CC(=O)Nc1ccc(O)cc1",
        "name": "Paracetamol",
        "terminal_cs_threshold": 1.5,
    },
)

ctx = cmd.execute("next")

cmd.execute(
    "try_precursors",
    {
        "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
        "reaction_type": "Schotten-Baumann acylation",
    },
)

cmd.execute(
    "commit",
    {
        "idx": 0,
        "reasoning": "Acylation with simple, accessible precursors.",
        "confidence": "high",
    },
)
```

这是一个协议层面的示例，而不是完整 benchmark 工作流。更多技术说明保存在 [usage notes](docs/usage-notes.md) 中。

<details>
<summary><strong>仓库结构</strong></summary>

- [main](main): 编排逻辑、会话逻辑、路线树、报告与命令接口
- [chem_tools](chem_tools): 具备化学约束的操作器与验证工具
- [tools](tools): 运行、分析、可视化及相关研究流程的辅助脚本
- [tests](tests): 当前验证与实验支持材料
- [plan](plan): 论文草稿、写作材料与论文准备资源

</details>

## 项目状态

- 活跃研究代码库
- 正在为面向 arXiv 的展示进行整理
- 文档正在持续清理中，但仓库仍处于日常活跃使用状态
- 核心工作流已经在使用中
- 尚未达到完全打磨后的开源发布状态
