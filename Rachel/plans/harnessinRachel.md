# Harness In Rachel

> 内部工程分析文档。
> 本文不是把 Rachel 重新命名为一个独立的 “harness 子工程”，而是从 harness 视角提炼 Rachel 现有系统里已经成型的编排、上下文、工具协议、化学门控与审计结构。

## 总论

Rachel 当初并不是以 “harness 工程” 这套命名来设计的，但从今天的代码形态回看，它实际上已经形成了一套非常完整的 harness-like 系统骨架。

如果把 harness 理解为“包裹模型的一层结构化运行环境”，它至少要回答五个问题：

1. 模型通过什么动作面与系统交互。
2. 当前状态由谁持有，怎样恢复。
3. 上下文怎样压缩、分页、逐步展开。
4. 候选动作怎样先试探、再提交，而不是直接污染主状态。
5. 领域约束怎样作为硬门控和软规则进入决策过程。

Rachel 对这五个问题都给出了系统级答案：

- `RetroCmd` 提供稳定的 JSON 命令协议。
- `RetroSession` 提供会话快照、沙盒保存和跨进程恢复。
- `RetrosynthesisOrchestrator` 提供队列调度、当前上下文、终止判定与提交流程。
- `ProposalContext` 提供 `compact / full / status / tree` 的分层上下文与窗口化摘要。
- `try_* -> select -> commit` 提供明确的沙盒纪律。
- `forward_validate`、`site_audit`、`fg_warnings`、`smart_cap` 与 `skill.md`/`experience.md` 一起构成领域约束层。

因此，更准确的说法不是“Rachel 里有一个显式名为 harness 的子目录”，而是“Rachel 本身已经内生出一个面向多步逆合成的领域化 harness”。它和通用 agent harness 的目标不同，但在编排、上下文规划、工具暴露、渐进披露、失败保留和可恢复性上，已经与 harness 的核心精神一脉相承，殊途同归。

## 阅读立场

本文坚持代码事实优先，并明确区分三类信息来源：

- 运行时硬实现：由 Python 代码直接执行和强制。
- 操作手册与经验先验：由 `skill.md`、`experience.md` 提供，约束 LLM/操作者的工作方式，但不一定全部被代码硬编码。
- 测试证据：由 `main/_test_sandbox.py`、`tests/test_compact_bond_identifiers.py`、`tests/test_site_anchor_audit.py` 等测试体现，说明哪些系统行为已经被显式验证。

这一区分很重要。Rachel 的强项恰恰在于它不是单靠“代码门控”或单靠“提示词规则”成立，而是把二者叠加成一个更完整的领域 harness。

## 1. Harness 视角下的 Rachel 定位

### 1.1 Rachel 不是一次性 route generator

Rachel 的根本差异，不在于它能不能输出一条路线，而在于它把路线规划从一次性文本生成改写成了一个持久化的决策过程。

根层文档已经把这一点说得很清楚：

- `README.zh-CN.md` 把 Rachel 描述为 `state -> action -> validation -> commit` 过程。
- `refs.md` 把项目定位为 “通过 JSON 命令接口让 LLM 与化学工具层交互，将目标分子逐步拆解为简单原料，并把决策过程保留为可回放、可审计的会话对象”。
- `workflow.md` 明确把系统拆成初始化、规划循环、终止、报告四个阶段。

从 harness 视角看，这意味着 Rachel 的主要职责不是“替模型算出答案”，而是“为模型提供一个可控、可恢复、可审计的运行框架”，让模型在框架内做局部决策。

### 1.2 更准确的定位

如果要给 Rachel 一个 harness 取向的工程定义，比较合适的说法是：

> Rachel 是一个面向多步逆合成的状态化决策 harness。
> 它把模型放在策略层，把化学工具放在动作与验证层，把会话、树结构、审计和导出放在编排层。

这个定义比“LLM 辅助逆合成程序”更准确，因为它强调的是：

- 模型不是唯一主角。
- 化学工具不是被动插件。
- 会话状态和失败轨迹不是附属日志，而是主系统的一部分。

### 1.3 Harness 对应物

| Harness 抽象 | Rachel 对应实现 | 作用 |
| --- | --- | --- |
| runtime shell | `RetroCmd` | 统一命令入口 |
| session container | `RetroSession` | 保存当前工作面、沙盒、队列、树、审计 |
| scheduler | `RetrosynthesisOrchestrator` | 负责 `next / prepare / commit / accept / skip / finalize` |
| context policy | `ProposalContext` | 决定默认给模型看什么、何时展开细节 |
| tool layer | `chem_tools/*` + `tools/llm_retro_platform.py` | 生成候选、化学分析、验证与提示 |
| durable memory | `RetrosynthesisTree` + `SynthesisAuditState` | 保留路线、决策、失败、保护策略 |
| artifact surface | `report / export / visualization` | 输出可检查结果 |

## 2. 编排与状态机

### 2.1 责任链

Rachel 的主编排链条非常清晰：

`RetroCmd -> RetroSession -> RetrosynthesisOrchestrator -> RetrosynthesisTree / SynthesisAuditState`

这条链对应了一个典型 harness 的四层分工。

#### `RetroCmd`

`main/retro_cmd.py` 是最外层命令面板。它把系统暴露成一组显式命令，而不是要求模型直接操作内部对象。这里的核心不是“封装得漂亮”，而是把动作边界固定下来。

#### `RetroSession`

`main/retro_session.py` 负责把“当前工作面”持久化成一个单文件会话对象。它持有：

- 当前目标与配置
- 当前活跃上下文
- 沙盒尝试历史与选中项
- 队列与 seen 集
- 合成树
- 审计状态

这使 Rachel 的运行单位不是某次提示词对话，而是一个完整的 session。

#### `RetrosynthesisOrchestrator`

`main/retro_orchestrator.py` 是真正的状态机内核。它维护：

- `_queue`
- `_seen`
- `_current_context`
- `_steps_executed`
- `tree`
- `audit_state`

同时负责：

- 从队列取下一个待处理分子
- 决定是否 quick-pass terminal
- 构建标准决策上下文
- 执行 commit
- 推进树与队列

#### `RetrosynthesisTree` 与 `SynthesisAuditState`

`main/retro_tree.py` 和 `main/retro_state.py` 则分别承担 durable memory：

- 前者保留“路线已经变成了什么”
- 后者保留“系统是怎么走到这里的”

这两个对象组合在一起，才让 Rachel 具备真正意义上的可回放性。

### 2.2 状态机长什么样

如果用最紧凑的工程视角重写 Rachel 主循环，它大致是：

```text
init
  -> next / prepare_next
  -> 生成 ProposalContext
  -> explore / try_*
  -> select
  -> commit 或 accept 或 skip
  -> 更新 queue / tree / audit
  -> 重复直到 finalize
```

这和普通脚本式 “输入目标 -> 输出路线” 的差别在于，Rachel 每一步都要求显式的状态转移。

### 2.3 关键状态行为

`RetrosynthesisOrchestrator.prepare_next()` 体现了几个很重要的 harness 思路：

- 若上一个 `_current_context` 未完成，会先放回队列，而不是静默丢失。
- quick-pass terminal 与 standard decision 是两条不同的 decision tier。
- `build_decision_context(smiles)` 只在需要标准决策时生成。
- 构建新上下文时，会把 `queue_preview`、`audit_state_summary`、`failed_attempts_for_current` 一起带进来。

这意味着 Rachel 的“当前一步”并不是孤立的，而是带着局部历史、全局进度和待办预览进入下一轮。

### 2.4 为什么这很像 harness

通用 harness 往往会把模型放进一个 “step loop” 中，要求模型一次只对当前任务面做一小步决策。Rachel 已经完全具备这一点，而且比很多通用 harness 更严格：

- 它的 step 不是抽象推理步，而是化学规划步。
- 它的状态不是自由文本，而是树、队列、审计和上下文对象。
- 它的调度不是聊天顺序，而是明确的 BFS 式展开。

### 2.5 测试证据

- `main/_test_sandbox.py` 通过 `prepare_next -> explore_bond -> try_disconnection -> commit_decision` 的完整流程，验证了 Rachel 的主决策循环不是概念图，而是可执行协议。
- `main/_test_session.py` 验证了 session 保存、加载、继续运行和 tree/report 导出之间的一致性。

## 3. 上下文规划与渐进式披露

### 3.1 Rachel 最像 harness context policy 的地方

如果只选 Rachel 中最有 “harness 味道” 的一部分，那大概率就是上下文规划。

Rachel 不是把所有分子信息、所有模板、所有前体、所有审计数据一股脑喂给模型，而是明确设计了上下文分层与按需展开机制。

这个机制主要由两部分组成：

- `tools/llm_retro_platform.py` 负责构建 full decision context
- `main/retro_orchestrator.py` 里的 `ProposalContext` 负责把 full context 压缩成多种视图

### 3.2 Full context 与 compact context 的边界

`build_decision_context(smiles)` 会收集：

- 分子基础信息
- 官能团摘要
- 复杂度分类
- `disconnectable_bonds`
- `fgi_options`
- `warnings`

但这只是 tool-layer full context，不等于默认给 LLM 的上下文。

`ProposalContext._build_compact_payload()` 会进一步做几件事：

- 只保留 `molecule / functional_groups / complexity / warnings` 的摘要信息
- 对 bond 候选只给 `bond_summary`
- 对 FGI 候选只给 `fgi_summary`
- 使用显式 `bond_offset / bond_limit / fgi_offset / fgi_limit`
- 保留 `bonds_omitted / fgi_omitted`
- 给出下一步动作提示

这里的关键不只是“压缩”，而是“保留索引稳定性并支持分页”。这已经非常接近成熟 harness 里的 context window policy。

### 3.3 渐进式披露的三个层级

Rachel 的渐进披露大体可以分成三层：

#### 第一层：compact summary

默认 `next()` 和 `context(detail="compact")` 只返回第一页窗口。这是默认工作面。

它让模型先看：

- 当前分子是什么
- 有哪些官能团
- 有多少 bond/FGI 候选
- 当前窗口里最值得看的候选是谁

但不直接暴露所有完整前体方案。

#### 第二层：local deep dive

当模型需要更深信息时，再用：

- `explore_bond(idx)`
- `explore_fgi()`

这一步才展开某个 bond 或某批 FGI 的完整前体信息。

值得注意的是，`explore_bond()` 还会额外返回 `bond_fg_context`，也就是围绕某个具体键位的局部官能团微审计。这是一种非常典型的 “局部放大镜” 设计。

#### 第三层：global status/tree view

当需要看全局时，再使用：

- `context(detail="status")`
- `context(detail="tree")`
- `tree`

这类视图不是用于日常每一步决策，而是用于状态确认、调试和全局复盘。

### 3.4 信息压缩不是削弱，而是去噪

`tools/llm_retro_platform.py` 和 `main/retro_orchestrator.py` 里都有一类非常重要的注释：保留旧逻辑作为注释，而当前逻辑改为更紧凑、更 LLM-friendly 的摘要。

这说明 Rachel 的上下文设计已经经历过一次工程上的认识升级：

- 不是信息越多越好
- 不是 raw SMARTS hit map 直接暴露给 LLM 就更专业
- 真正重要的是给模型“足够支撑决策”的最小结构化信息

这种设计很符合 harness 的成熟形态：让上下文成为被治理的资源，而不是自然泄漏的内部实现。

### 3.5 `actual_bond_idx` 与 `role_pair`

`tests/test_compact_bond_identifiers.py` 说明 Rachel 不只是压缩信息，还在 compact 层保留了最关键的辨识锚点：

- `actual_bond_idx`
- `role_pair`

这两个字段的意义很大：

- `actual_bond_idx` 把 abstract bond proposal 对回真实 RDKit bond
- `role_pair` 把原子环境压缩成模型更容易理解的化学角色

这说明 Rachel 的 context policy 不是“少给”，而是“少给噪声，多给辨识度高的锚点”。

### 3.6 Harness 对应物

| Harness 概念 | Rachel 中的对应物 |
| --- | --- |
| context policy | `ProposalContext` |
| default context window | `compact` |
| paged retrieval | `bond_offset / bond_limit / fgi_offset / fgi_limit` |
| local drill-down | `explore_bond / explore_fgi` |
| global debugging view | `status / tree / full` |

### 3.7 测试证据

- `main/_test_sandbox.py::test_layered_context` 直接验证 `compact` 与 `full` 的差异。
- `tests/test_compact_bond_identifiers.py::test_compact_context_supports_windowed_bond_and_fgi_summaries` 验证 compact window 的分页契约。
- `tests/test_compact_bond_identifiers.py::test_explore_bond_adds_local_bond_fg_context_without_reinjecting_full_map` 验证局部展开不会把旧式完整噪声重新注入默认上下文。

## 4. 工具暴露与动作协议

### 4.1 Rachel 的工具面不是隐式方法，而是显式协议

从 harness 视角看，Rachel 的另一条核心主线是它把系统能力暴露成了一个非常明确的 action surface。

`main/retro_cmd.py` 中的 `RetroCmd.execute()` 通过 dispatch 暴露出一整组标准动作：

- `init`
- `next`
- `context`
- `explore`
- `explore_fgi`
- `try_bond`
- `try_fgi`
- `try_precursors`
- `sandbox_list`
- `sandbox_clear`
- `select`
- `commit`
- `accept`
- `skip`
- `tree`
- `status`
- `finalize`
- `report`
- `export`
- `smart_cap`
- `custom_cap`

这已经不是“提供几个函数调用”，而是一套完整的协议面。

### 4.2 JSON-in / JSON-out 的意义

`retro_cmd.py` 顶部注释明确把自己定义为 “标准化的 JSON-in / JSON-out 命令层”。这件事很关键。

它意味着 Rachel 的外部模型接口遵循的是一种工具协议，而不是：

- 直接 import 内部对象
- 依赖 REPL 状态
- 把半结构化文本当作唯一交互渠道

对于 harness 来说，这种协议有三个直接好处：

- 动作边界稳定
- 结果格式稳定
- 中间状态可以被文件系统和外部调度器接管

### 4.3 命令协议的分层

这些命令并不是平铺的，它们实际上形成了五类动作：

#### 会话管理

- `init`
- `status`
- `finalize`

#### 决策工作面

- `next`
- `context`
- `tree`

#### 局部展开

- `explore`
- `explore_fgi`

#### 沙盒试探与选择

- `try_bond`
- `try_fgi`
- `try_precursors`
- `sandbox_list`
- `sandbox_clear`
- `select`

#### 提交与导出

- `commit`
- `accept`
- `skip`
- `report`
- `export`

这种分层让 Rachel 的动作面具有良好的工程可解释性：模型不是直接“做路线”，而是在一套已定义好的动作语法里推进路线。

### 4.4 `cmd.json` 模式的意义

`retro_cmd.py` 还保留了 `cmd.json -> result.json` 的文件桥接模式。

这一点说明 Rachel 的作者当时已经在做一件非常 harness 的事：把模型与运行时解耦。

也就是说，即便不在同一个 Python 进程里，模型依然可以通过：

- 写命令文件
- 触发执行
- 读取结果文件

来驱动整个系统。这是非常典型的 protocol-oriented harness 设计。

### 4.5 为什么这比“工具很多”更重要

Rachel 的价值不只是工具多，而是这些工具被组织成了一个有顺序、有前后约束的协议。

例如：

- 不能在没有活跃上下文时直接 `commit`
- `try_*` 与 `commit` 被故意拆开
- `select` 明确指向某次沙盒尝试
- `context` 与 `explore` 负责信息获取，不负责状态污染

这说明 Rachel 的工具暴露不是功能拼盘，而是行为编排的一部分。

### 4.6 Harness 对应物

| Harness 概念 | Rachel 中的对应物 |
| --- | --- |
| tool contract | `RetroCmd.execute()` dispatch |
| action schema | 每个 `cmd_*` 的参数与返回结构 |
| external adapter | `cmd.json / result.json` 模式 |
| controlled side effects | `try_*` 读写沙盒，`commit` 才写主树 |

## 5. 沙盒纪律与提交门控

### 5.1 Rachel 最重要的工程纪律之一

Rachel 最像成熟 harness 的另一个地方，是它对 “试探” 和 “提交” 之间边界的强调。

`refs.md` 已经把这点概括为：

> `try_*` 不写入主树，`commit` 后才持久化。

这不仅是实现细节，而是 Rachel 整个规划哲学的核心纪律。

### 5.2 沙盒在 Rachel 中是什么

`RetroSession` 持有：

- `_sandbox_attempts`
- `_sandbox_selected`

任何一次：

- `try_bond`
- `try_fgi`
- `try_precursors`

都会生成一个 `SandboxResult`，被追加到 session 级沙盒历史里。

从 harness 角度看，这个沙盒就是 Rachel 的 ephemeral workspace：

- 可反复比较
- 可擦除
- 与主树隔离
- 可在 gate fail 后保留

### 5.3 提交前的显式选择

Rachel 不是默认拿最后一次尝试去提交，而是要求：

1. 先做多次 `try_*`
2. 用 `sandbox_list` 看历史
3. 用 `select` 显式选中方案
4. 再 `commit`

这一点很重要，因为它把“我以为要提交的是哪个方案”变成了结构化状态，而不是依赖上下文记忆。

### 5.4 commit 不是简单写树

`RetroSession.commit()` 在真正调用 `orch.commit_decision()` 之前，还会做一层额外 gate：

- `_evaluate_site_retention_gate()`

这个 gate 的作用是：

- 对 `llm_proposed` 的 same-core 手写前体做位点保真审计
- 如果同骨架 proposal 缺少完整 site audit 信息，则阻断提交
- 如果出现 anchor drift 或 changed sites 过多，则阻断提交
- gate fail 后记录失败，但不销毁活跃沙盒

这说明 Rachel 的提交门控不是一个单点函数，而是：

`Sandbox attempt -> session-level gate -> orchestrator commit -> tree/audit update`

### 5.5 gate fail 后保留什么

`RetroSession._record_gate_failed_attempt()` 会把失败写入：

- `audit_state.record_failure`
- `audit_state.record_decision(... outcome="gate_failed")`

而且会保留当前沙盒状态。

这是非常重要的 harness 特征：失败不是单纯报错退出，而是变成系统记忆的一部分。

### 5.6 自动附带保护建议

`RetroSession.sandbox_try()` 里还有一个很典型的领域化 harness 细节：

如果 `forward_validation.hard_fail_reasons` 中出现 `forbidden_fg`，系统会尝试调用 `suggest_protection_needs()`，把可执行的保护建议附在 attempt 上。

这代表 Rachel 的沙盒并不只是 “判你不行”，而是开始往 “为什么不行、下一步怎么办” 进化。

### 5.7 prepare_next 为什么会清空沙盒

`RetroSession.prepare_next()` 在取新分子前会清空当前沙盒。

这看似普通，实际也是一种纪律设计：

- 沙盒属于当前活跃 decision context
- 一旦切换 decision target，旧沙盒默认失效
- 避免不同分子的候选历史污染当前选择

也就是说，Rachel 的沙盒是 session 级对象，但作用域是当前 decision surface。

### 5.8 Harness 对应物

| Harness 概念 | Rachel 中的对应物 |
| --- | --- |
| ephemeral workspace | `_sandbox_attempts` |
| proposal comparison set | `sandbox_list` |
| explicit commit target | `sandbox_select` |
| safety gate before state mutation | `_evaluate_site_retention_gate` |
| failed attempt retention | `audit_state.record_failure` |

### 5.9 测试证据

- `main/_test_sandbox.py` 验证 `explore + sandbox + commit` 的工作流。
- `tests/test_site_anchor_audit.py::test_retro_session_commit_blocks_bad_same_core_precursor_and_keeps_sandbox` 验证 gate fail 时沙盒被保留，审计状态被更新。
- `tests/test_site_anchor_audit.py::test_retro_session_commit_allows_same_core_single_site_change_even_when_scaffold_alignment_is_below_one` 验证 Rachel 的 gate 不是粗暴依赖全局 scaffold score，而是更细粒度地使用 site audit。

## 6. 化学专业性如何进入 harness

### 6.1 Rachel 不是通用 agent harness

Rachel 最有价值的一点，不是它拥有一套通用模型编排骨架，而是这套骨架从一开始就被化学专业性塑形。

它不是 “先有通用 harness，再把化学工具挂上去”，而更像：

- 编排层从一开始就知道自己在组织的是 retrosynthesis step
- 上下文层从一开始就在压缩化学关键信息
- 门控层从一开始就服务于化学可行性与位点保真

因此，Rachel 更准确的工程类别是：

> 嵌入化学审计规则的领域 harness

### 6.2 代码中被硬实现的专业层

#### `forward_validate.py`

`chem_tools/forward_validate.py` 提供的是 Rachel 最核心的硬门控之一。

`validate_forward()` 不是一个模糊打分器，而是顺序执行的验证栈：

1. atom balance
2. template forward execution
3. scaffold alignment
4. bond topology
5. reaction-specific plausibility
6. functional-group compatibility

并把若干条件定义为 hard gate，例如：

- `severe_imbalance`
- `skeleton_imbalance`
- `scaffold_not_aligned`
- `bond_topology_violation`
- `reaction_specific_violation`
- `forbidden_fg`

这说明 Rachel 的“可提交性”不是由模型主观判断决定，而是由化学验证栈显式介入。

#### `site_audit.py`

`chem_tools/site_audit.py` 进一步把普通 scaffold alignment 无法解决的问题单独拉了出来：同骨架 proposal 的位点保真。

它会：

- 识别 same-core major precursor
- 建立 scaffold atom mapping
- 比较 changed sites 与 preserved sites
- 生成 summary、site_rows、changed_site_count

这正是 Rachel 在复杂杂环和药化样分子上相对“纯模板系统”更专业的一点。

#### `smart_cap.py`

`chem_tools/smart_cap.py` 则把 bond environment 转成更像真实前体的 capping 建议。

它覆盖的不是抽象模板名，而是现实中常见的合成连接方式，例如：

- `Suzuki`
- `amide bond formation`
- `Buchwald-Hartwig`
- `SNAr / Ullmann`
- `Grignard`

这使 Rachel 的候选空间不只来自模板扫描，也来自基于键环境的推断。

#### `fg_warnings.py`

`chem_tools/fg_warnings.py` 则提供：

- 官能团冲突检测
- 反应兼容性检查
- 保护需求建议
- 脱保护安全性检查

这让 Rachel 的 guardrail 不只是“这个 proposal 合不合理”，还包括“这一步与现有功能团生态是否冲突”。

### 6.3 规则文档中的专业层

Rachel 的专业性并不只存在于代码里，还系统性地写进了 `skill.md` 和 `experience.md`。

#### `skill.md` 的作用

`skill.md` 更像 Rachel 的操作宪法。它规定了：

- 决策优先级
- commit 前必须审计什么
- 什么时候允许手写前体
- 什么情况下必须显式写保护/脱保护
- 为什么 `pass=true` 只能作参考
- 为什么要优先收敛设计
- 什么时候可以接受 advanced terminal

这里的重点是：这些规则不是普通 prompt 润色，而是在给模型规定 “怎样在这个 harness 里行动才算合格”。

#### `experience.md` 的作用

`experience.md` 则更像经验层排序器。它不一定都是硬规则，但会持续影响候选排序与怀疑方式，例如：

- 先数骨架原子，再谈反应类别
- 成熟药化杂环优先保骨架
- 复杂片段连接优先考虑收敛路线
- 高反应性手柄默认晚装
- 每一步都要追关键碳来源

这相当于在硬门控之外，再给 harness 增加一层资深化学家式的 heuristic policy。

### 6.4 代码门控与规则手册如何协同

Rachel 的专业性不是单一来源，而是三层叠加：

| 层级 | 载体 | 作用 |
| --- | --- | --- |
| tool-layer chemistry | `chem_tools/*` | 生成候选、验证可行性、检测冲突 |
| orchestration gates | `RetroSession` / `RetrosynthesisOrchestrator` | 决定何时允许写入主树 |
| operator policy | `skill.md` / `experience.md` | 规定模型如何解释、排序、拒绝候选 |

这使 Rachel 的化学专业性不是“附会在输出文本中的行话”，而是已经进入 harness 结构。

### 6.5 需要明确区分的边界

这里必须强调一个工程判断：

- `forward_validate.py`、`site_audit.py`、`fg_warnings.py`、`smart_cap.py` 的行为是运行时硬实现。
- `skill.md`、`experience.md` 中的很多优先级和偏好，是操作规则，不等于所有条目都被代码硬编码。

也就是说，Rachel 的专业性很强，但它的一部分是 code-enforced，一部分是 policy-enforced。二者共同构成系统，而不是互相替代。

### 6.6 Harness 对应物

| Harness 概念 | Rachel 中的对应物 |
| --- | --- |
| domain guardrails | `skill.md` + `experience.md` |
| hard validators | `forward_validate.py` |
| site-specific safety gate | `site_audit.py` + session gate |
| tool-augmented proposal shaping | `smart_cap.py` |
| reaction compatibility layer | `fg_warnings.py` |

## 7. 可审计性、可恢复性与外显产物

### 7.1 会话恢复不是附带功能，而是主设计

`RetroSession.to_dict()` 导出的内容非常完整：

- `target`
- `config`
- `status`
- `current`
- `queue`
- `seen_smiles`
- `tree`
- `audit_state`

其中 `current` 还复用了 live `compact` 生成逻辑，并附带当前沙盒。

这说明 Rachel 不只是把“结果”持久化，而是把“当前工作面”也持久化。对于 harness 来说，这是非常重要的成熟特征。

### 7.2 跨进程恢复当前 decision surface

`RetroSession.load()` 不只是把树和状态读回来，它还会：

- 恢复 `_sandbox_attempts`
- 恢复 `_sandbox_selected`
- 重建 `_current_context`
- 重新生成 `decision_context`

也就是说，一个中断的 decision surface 可以被继续操作，而不是只能从树结构里推断回去。

这让 Rachel 的会话恢复是“继续工作”，不是“从头复盘”。

### 7.3 Tree 记录“路线是什么”

`RetrosynthesisTree` 使用双类型节点图：

- `MoleculeNode`
- `ReactionNode`

并在反应节点里保存：

- `reaction_smiles`
- `reaction_type`
- `template_evidence`
- `llm_decision`
- `forward_validation`

这意味着主树不仅是拓扑结构，也保留了提交时的证据和解释。

### 7.4 Audit 记录“为什么变成这样”

`SynthesisAuditState` 则记录：

- `decision_history`
- `failed_attempts`
- `protections`
- strategic summary

从工程意义上说：

- tree 是 result memory
- audit_state 是 process memory

Rachel 同时保留了二者，因此能同时回答：

- “现在路线长什么样？”
- “这条路线为什么会变成现在这样？”

### 7.5 导出产物是 artifact pipeline

`retro_report.py`、`retro_visualizer.py`、`retro_output.py` 构成了一条完整的 artifact pipeline。

`export_results()` 会导出：

- `SYNTHESIS_REPORT.html`
- `SYNTHESIS_REPORT.md`
- `report.txt`
- `tree.json`
- `tree.txt`
- `terminals.json`
- `visualization.json`
- `session.json`
- `images/`

这不是简单的“出几份报告”，而是在把决策过程变成多个面向不同用途的可视化与数据产物：

- 给人读的文本报告
- 给前端和分析脚本用的 JSON
- 给复盘用的 session snapshot
- 给展示用的 HTML/图像

从 harness 角度看，这是一种很强的 artifact externalization 能力。

### 7.6 为什么 artifact 很重要

许多系统只保留最终答案，失败过程和中间状态在运行结束后就消失了。

Rachel 相反：

- 保留 failed attempts
- 保留 rejected alternatives
- 保留 selection reasoning
- 保留 forward validation
- 保留 session snapshot

这使它天然适合：

- 方法分析
- 错误归因
- benchmark 复盘
- prompt/protocol 迭代

### 7.7 Harness 对应物

| Harness 概念 | Rachel 中的对应物 |
| --- | --- |
| durable memory | `RetrosynthesisTree` + `SynthesisAuditState` |
| resumable execution | `RetroSession.save/load` |
| run artifact bundle | `export_results()` |
| replayable trace | `session.json` + `decision_history` + `tree.json` |

## 关键接口与契约摘要

这一节把本文前面分散提到的关键接口压缩成一页工程摘要。

### `RetroCmd` 命令协议

`RetroCmd` 是 Rachel 的标准 action surface。它的核心契约是：

- 输入为 `command + args`
- 输出为结构化 dict
- 未知命令显式报错
- 所有动作通过固定 dispatch 进入 session

这使模型和外部系统不需要知道内部类的组织方式，只需要知道协议。

### `RetroSession` 会话与快照契约

`RetroSession` 的核心契约是：

- 当前工作面必须可保存
- session 恢复后应尽可能回到原 decision surface
- sandbox 与 tree 同属 session，但职责不同
- `current` 的 compact 输出应与 live context 保持一致，不允许 session snapshot 与实时上下文漂移

### `ProposalContext` 分层上下文契约

`ProposalContext` 的核心契约是：

- `compact` 是默认工作面
- `full` 是重型视图，不应作为默认
- `status` / `tree` 用于全局状态确认
- compact 必须支持显式窗口参数
- local deep dive 通过 `explore_bond / explore_fgi` 完成

### `RetrosynthesisTree` 与 `SynthesisAuditState` 的记忆契约

这两个对象分别保证：

- 主树记录 accepted route structure
- 审计状态记录 decisions, failures, protections
- 二者一起组成 durable memory，而不是互相替代

### 导出语义契约

Rachel 导出的几个关键 JSON 各自承担不同语义：

- `session.json`: 可恢复工作态
- `tree.json`: 已提交主树结构
- `visualization.json`: 前端/图结构消费面
- `terminals.json`: 起始原料清单

这说明 Rachel 的 artifact 并不是同一份数据的重命名，而是面向不同消费场景的结构化导出。

## Rachel 的 Harness 优点

- 它把模型运行从一次性回答变成了状态化、多轮、可恢复的规划过程。
- 它的动作面稳定且明确，模型与系统的边界清楚，不依赖隐式上下文猜测。
- 它对上下文做了真正的治理，默认摘要、窗口化分页、局部展开、全局视图各有边界。
- 它把试探与提交严格分开，让候选探索不会直接污染主树。
- 它把化学专业性深嵌到 harness 里，而不是只把化学术语写进输出文本。
- 它保留失败、选择理由、验证结果和导出产物，天然适合复盘、分析和迭代。
- 它同时拥有 result memory 与 process memory，这比只保存最后一棵树更有研究价值。

## Rachel 还没有显式产品化的 Harness 部分

Rachel 虽然已经具备很强的 harness 骨架，但它目前仍更像“内生于研究系统的 harness”，而不是显式产品化的通用 harness 框架。

主要体现在以下几点：

- 仓库里没有独立命名的 `harness/` 模块或单独包装的公共子包。
- 编排、规则、实验脚本、论文材料和数据资产仍共处于同一个研究仓库里。
- 很多 harness 语义是通过 `retro_*`、`session`、`tree`、`audit` 这些领域命名来表达，而不是以通用 runtime 抽象显式抽离。
- 上下文策略和动作协议已经比较稳定，但还没有被声明为独立版本化的公共 API。
- 领域规则分散在代码、`skill.md`、`experience.md`、测试和文档中，尚未统一抽象成单一的 policy subsystem。
- 它的设计明显服务于 retrosynthesis 这一特定任务，而不是追求跨任务通用性。

但这并不是缺点。恰恰相反，这正说明 Rachel 的强项不在“抽象得像个平台”，而在“把特定领域的问题组织成一套可控、可验证、可恢复的运行机制”。

也可以说，Rachel 的 harness 不是被产品经理先定义出来的，而是在长期解决真实路线规划问题的过程中自然长出来的。

## 证据锚点

下面这些文件是本文结论最核心的证据来源：

- `main/retro_cmd.py`
  - 命令协议、JSON-in / JSON-out、`cmd.json` 桥接模式
- `main/retro_session.py`
  - session 快照、沙盒管理、gate fail 保留、save/load 恢复
- `main/retro_orchestrator.py`
  - `ProposalContext`、BFS 调度、`prepare_next`、`commit_decision`
- `main/retro_tree.py`
  - `MoleculeNode / ReactionNode`、TemplateEvidence、LLMDecision
- `main/retro_state.py`
  - decision/failure/protection 的 durable audit memory
- `tools/llm_retro_platform.py`
  - `build_decision_context`、FG summary 压缩、bond/FGI 上下文生成
- `chem_tools/forward_validate.py`
  - 正向验证栈与 hard fail reasons
- `chem_tools/site_audit.py`
  - same-core site retention 审计
- `chem_tools/smart_cap.py`
  - bond environment 到 precursor style proposal 的映射
- `chem_tools/fg_warnings.py`
  - 功能团冲突、兼容性与保护建议
- `skill.md`
  - Rachel 的操作宪法与审计要求
- `experience.md`
  - 经验先验、常见误判与排序偏好

与本文最相关的测试包括：

- `main/_test_sandbox.py`
  - layered context、explore、sandbox、commit 的主流程验证
- `tests/test_compact_bond_identifiers.py`
  - compact context 的锚点信息、分页窗口与局部放大镜验证
- `tests/test_site_anchor_audit.py`
  - same-core site gate、gate fail 保留沙盒、commit 阶段门控验证

## 最终判断

如果从名称上看，Rachel 里没有一个显式叫 harness 的工程。

但如果从系统结构看，Rachel 已经具备了一个成熟领域 harness 的关键构件：

- 明确动作协议
- 持久化会话状态
- 分层上下文策略
- 沙盒试探与门控提交
- 领域硬验证
- 失败保留与审计记忆
- 外显 artifact 导出

所以，Rachel 更合适的理解方式不是“一个带很多工具的逆合成脚本集合”，而是“一个已经长成 harness 形态的多步逆合成运行框架”。

只是这个 harness 目前仍以内生方式存在于 Rachel 之中，而没有被单独产品化命名出来。
