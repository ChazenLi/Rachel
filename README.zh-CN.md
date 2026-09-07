[English](README.md) | [简体中文](README.zh-CN.md)

# Rachel

**由 LLM 主导的多步逆合成：从连续化学决策到完整路线。**

[工作原理](#how-it-works) · [研究结果](#research) · [完整路线图谱](#route-atlas) · [开始使用](#get-started) · [许可证](#license)

每一次逆合成选择，都会改变接下来需要制备的分子。一次局部合理的断开，可能留下更难合成的前体，
引入新的选择性要求，甚至需要重新考虑整条路线的策略。要完成多步规划，化学判断就必须随着问题的变化持续进行。

**Rachel 为通用大语言模型提供了一个能够持续开展这项工作的化学规划环境。**
模型提出转化、比较备选方案，并在必要时修订路线策略；Rachel 提供分子上下文，
在计算沙箱中执行候选转化、返回验证证据，并将已接受的步骤带入下一轮决策。
逐步形成的路线及其背后的选择，都保留为可检查的记录。

本仓库包含 **Rachel-beta 运行时**，以及论文采用的 **PaRoutes120 与 RF25 完整路线图谱（Route Atlas）**。
你可以直接浏览已有结果，也可以将运行时接入能够调用工具的 LLM 智能体，开展路线构建与复审。

<a id="how-it-works"></a>
## 工作原理

![Rachel 系统总览：LLM 维持路线策略，分子上下文、候选执行和验证证据共同支持连续决策。](docs/assets/rachel-workflow.png)

*论文中的系统总览。已接受的转化更新路线状态，也改变下一轮决策面对的分子问题。
上半部分对照不同的路线组织方式，不表示性能排名。*

LLM 负责合成策略，环境则将每个选择的后果明确保留下来。处理当前分子时，
模型能够同时参考整体计划、结构约束和已经作出的决策。

1. **提出下一步转化。** 模型可以考察预定义的断键与官能团转化操作，也可以提出自己的前体组合与反应方案。
2. **考察化学证据。** 沙箱生成候选前体，并返回原子来源、反应位点、拓扑及兼容性等检查结果。
   模型据此比较候选、处理警告或修订提议；硬性阻断的候选不能提交。
3. **从接受的状态继续。** 提交转化后，路线树随之更新。仍需合成的前体成为后续规划问题，
   整体策略也可以随着路线的展开而修订。

化学工具为决策提供依据，LLM 或化学家负责解释这些依据。
检查通过说明验证器在其能力范围内取得了支持，尚缺少的证据也会保留在规划记录中。

<a id="research"></a>
## 研究结果

我们的研究考察：当前序选择不断改变后续分子问题时，通用 LLM 能否持续组织整条路线的化学决策。
GPT-5.5 在 Rachel 中从目标结构出发开展规划，不接收参考路线或路线级解答。

| 数据集 | 设置 | Rachel 严格闭合 |
| --- | --- | ---: |
| PaRoutes120 | 从 PaRoutes 中选取的 120 个目标 | 111 / 120 |
| RF25 | 单独分析的 25 个困难目标，来自分子探针与化学生物学研究 | 24 / 25 |

**严格闭合**要求路线构建完成，且每个终端前体均在规划结束后通过独立来源解析。
这一来源审计的信息没有提供给规划者。这些数字描述计算规划结果，实验室合成仍待检验。

研究进一步考察了路线如何形成，以及连接这些终端与目标的化学：

- **路线中的化学质量。** 在四种方法均闭合的 58 个共同 PaRoutes 目标上，
  前向模型的反应级评估与方法身份盲化的整路线评价提供了互补证据。
  两个整路线评估器均给予 Rachel 最高的平均总体评分，其中 Opus 5 对 Rachel 与 direct LLM 的评分接近。
  前向模型支持仍不完整，RF25 尤其如此。
- **随路线演化的决策。** 记录的轨迹显示，模型提出的化学持续进入路线的早期和后期位置；
  合成计划在过程中发生修订，后续步骤也能够对应到更新后的计划。
- **持续规划的作用。** 将 LLM 决策替换为固定策略后，局部化学仍能执行，
  但 PaRoutes 的严格闭合降至 6–15/120。限制规划支持功能同样降低了 RF25 的闭合数量。
  这些对比考察的是整个系统的行为。

公开路线记录保留了研究实际采用的版本，不能假定后续加入运行时的功能在所有记录中都已启用。

<a id="route-atlas"></a>
## 浏览完整路线图谱

Route Atlas 将队列结果对应到每一条具体分子路线。所有目标和方法位置均予以保留，
包括部分完成的路线、未取得输出的结果，以及 PaRoutes 的参考路线。

| 数据集 | 对比范围 | 离线页面 | JSON 记录 |
| --- | --- | --- | --- |
| PaRoutes120 | 120 个目标 × 8 种方法，另含 120 条 Reference 路线 | [PaRoutes120.html](data/route-atlas/PaRoutes120.html) | [PaRoutes120.json](data/route-atlas/data/PaRoutes120.json) |
| RF25 | 25 个目标 × 4 种方法 | [RF25.html](data/route-atlas/RF25.html) | [RF25.json](data/route-atlas/data/RF25.json) |

使用 **Code > Download ZIP** 下载并解压仓库，或通过 Git 克隆，
然后在浏览器中打开任一 HTML 文件。页面所需内容均已内置，
不需要服务器、网络连接、API key 或 Python 环境。GitHub 文件页提供下载，不直接运行交互页面。

选择一个目标后，可同时对照最多四种方法。分子结构和反应步骤关联已有的决策理由、
末端来源审计及匹配的计算评价。页面支持中英文界面、路线缩放和 JSON 导出。
PaRoutes120 默认并排展示 Rachel 与 Reference，RF25 默认展示 Rachel 与 Direct LLM。

数据包包含 **1,180 个对比位置**：1,060 个方法位置，以及单独列出的 120 个 Reference。
Reference 是数据集对照，不计入方法闭合数量。PaRoutes120 包含 Rachel、Direct LLM、
AOT*、SyntheLite、PaRoutesModel + Retro*、RootAligned + Retro*、
RootAligned + MCTS 和 LocalRetro + Retro*；RF25 包含前四种方法。

[数据说明](data/route-atlas/README.md) ·
[路线索引](data/route-atlas/ROUTE_INDEX.csv) ·
[字段定义](data/route-atlas/DATA_DICTIONARY.md) ·
[来源目录](data/route-atlas/SOURCE_CATALOG.csv)

<a id="get-started"></a>
## 开始使用

### 浏览已有结果

直接使用上方的 Route Atlas，无需模型访问权限或安装环境。

### 运行 Rachel

在已安装 Git 和 Conda 的终端中执行：

```bash
git clone https://github.com/ChazenLi/Rachel.git
cd Rachel
conda env create -f environment.yml
conda activate rachel-v2
```

所附环境指定 Python 3.11、RDKit 及配套依赖。公开运行时在 Windows 上准备，
尚未由本发行版验证其他平台。

让具备工具调用能力的 LLM 智能体访问本仓库，并读取
[Rachel 技能与命令协议](Rachel/skill.md)。模型接入由宿主智能体负责，
Rachel 通过 Python 命令接口提供化学工具和会话状态。

下面的代码创建独立会话，并查看初始规划上下文：

```python
from pathlib import Path
from tempfile import mkdtemp

from Rachel.tools.logged_runner import LoggedRetroCmd

run_dir = Path(mkdtemp(prefix="rachel_"))
cmd = LoggedRetroCmd(str(run_dir / "session.json"))

cmd.execute("init", {
    "target": "CC(=O)Nc1ccc(O)cc1",
    "name": "paracetamol",
})
print(cmd.execute("next", {}))
print(cmd.execute("context", {"detail": "compact"}))
print(f"Session and command logs: {run_dir}")
```

完整运行还需继续提出化学方案、进行沙箱检查、明确选择并更新路线，直至完成最终检查。
[命令协议](Rachel/skill.md)说明了这一循环及每次决策所需的证据。
`LoggedRetroCmd` 会在持久会话旁记录每条命令的输入和输出。

<a id="working-with-routes"></a>
## 路线复审与知识扩展

**检查与继续规划。** 会话状态、候选尝试和已提交反应分别保存。
`review_terminal` 可以重开末端，继续分解；`review_node` 可以从已经结束的路线建立独立变体。
因此，可以围绕某个分支开展复审，同时保留原路线。

**导出研究记录。** 报告和结构化输出包含会话、路线树、末端清单、分子可视化及知识配置记录，
既方便化学审阅，也可用于后续分析。

**扩展化学上下文。** 运行时支持基础、团队与项目层的版本化知识包，提供提示、经验和化学资源。
每个会话固定采用的知识配置，新增知识仍受运行时的验证与状态约束。
记录的实验事实可以用于形成尚未激活的知识草稿，经专家审阅并明确发布后，
新版本才能进入后续会话。

| 仓库入口 | 内容 |
| --- | --- |
| [Rachel/skill.md](Rachel/skill.md) | 智能体指令、公开命令与规划协议 |
| [Rachel/main](Rachel/main) | 会话编排、路线状态、报告与命令接口 |
| [Rachel/chem_tools](Rachel/chem_tools) | 分子分析、候选操作与验证 |
| [Rachel/knowledge](Rachel/knowledge) | 知识资源、版本化配置与知识包命令行工具 |
| [Rachel/tools](Rachel/tools) | 命令日志、末端审计与可视化工具 |
| [data/route-atlas](data/route-atlas/README.md) | 论文完整路线对比及对应证据记录 |

## 展望

逆合成为研究大模型如何接续化学决策提供了一个具体场景。更广泛的方向，
是将分子设计、合成规划与实验证据联系起来：一个候选分子难以制备，可以成为重新设计它的依据；
一项意外结果，也可能改变下一步的假设或实验。Rachel 使其中一部分交互能够在计算路线规划中被检查。
将其进一步扩展为可靠的实验研究协作，仍是未来的目标。

<a id="license"></a>
## 许可证与数据归属

Rachel 原创材料采用 **CC BY-NC 4.0** 许可。
非商业分享与改编须保留署名，商业用途不在该许可范围内。详见 [LICENSE](LICENSE)。

PaRoutes 参考数据来自
[PaRoutes，Zenodo 记录 7341155](https://doi.org/10.5281/zenodo.7341155)。
第三方材料保留其上游许可条件；页面内置库的声明见
[THIRD_PARTY_NOTICES.txt](data/route-atlas/THIRD_PARTY_NOTICES.txt)。
复用路线结果时，请保留数据集与来源标识，并记录所使用的仓库提交版本。
