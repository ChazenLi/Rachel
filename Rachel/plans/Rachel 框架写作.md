# Rachel 框架写作

## 1 定位：

形式化空间推理搜索----等intro来定一个主故事方向

## 2 实验资料准备：

### 1     单步数据库（USPTO50K+PAROUTES）来论证形式模板的化学空间有效效果与质量：

![USPTO 模板合理性分析](E:\Python\skills\Rachel\tests\baselines\results\collections_unified\01_uspto_btf_legacy_publication\figures\fig1_match_tier_waterfall.png)

上图是基于USPTO50K数据集进行的形式化模板效果验证，对于50K数据集中的单步反应rs>>ps，使用形式化模板来处理ps，得出可能的前体空间rs集合中包含ground truth，以及化学等价性的结果（其中化学等价性取非常严格>0.9，这部分是使用rdkit中的tanimoto进行相似度计算打分，可以作为参考）；

![](E:\Python\skills\Rachel\tests\baselines\results\collections_unified\01_uspto_btf_legacy_publication\figures\fig2_per_dataset_stacked.png)

上图是对整个USPTO50K数据集进行的形式化模板后处理的效果验证，可以看到整个目前形式化模板robust效果不错

同样对Paroutes中的n1\n5两万条多步合成数据全部进行单步拆分之后进行形式化模板质量验证分析：![](E:\Python\skills\Rachel\tests\baselines\results\collections_unified\03_paroutes_btf_bucket_reference_publication\figures\fig_pub_01_waterfall.png)

通过单步数据集实验来说明：

​	Rachel的形式化模板基底能够给出有效且高质量的化学空间；

​	Rachel的形式化模板可以初步实现受限空间约束搜索的效果；

​	仍然存在边界单相对数量很少可以通过完善形式化模板来丰富；

![](E:\Python\skills\Rachel\tests\baselines\results\collections_unified\05_conclusion\figures\conclusion_core3_composite_publication.png)

目前存在的问题：

​		如何和通用的模型进行效果的比较？整个地方确实不好进行对比，形式化模板对应的是一个潜在化学空间，不是明确的top-k形式；

​		后续这样写足够吗？从逻辑思考上是否足够？

----------



### 2     500个单步数据论证分析prompts的compact与full筛选实验

由于整个Rachel是面向LLM的化学形式化推理工作，故而每一个单步反应的prompts需要进行设计和优化，从而确保在效果、质量与上下文负载之间进行平衡取舍，下图是我们采用两种不同的反应prompts构建方式进行上下文工程在单步反应中的效果评价与质量构建得出的结果，full的prompts长度通常在compact的1.2x到2x之间；![](E:\Python\skills\Rachel\tests\E2\results\btf1to10_s500_per_question\llm_ready_blind_v2_split\comparison_compact_vs_full_20260307\fig02_tanimoto_bucket_comparison.png)

在从USPTO50K数据集中选取的500个单步数据集中，我们统计了llm的选择趋势：

![](E:\Python\skills\Rachel\tests\E2\results\btf1to10_s500_per_question\llm_ready_blind_v2_split\comparison_compact_vs_full_20260307\fig03_llm_choice_comparison.png)

此外对于目前的结果我们还需要确保单步输出的smiles合法性以及可以被解析的准确度：

![](E:\Python\skills\Rachel\tests\E2\results\btf1to10_s500_per_question\llm_ready_blind_v2_split\comparison_compact_vs_full_20260307\fig04_solved_exact_rate_comparison.png)

虽然compact在上下文管理之间更有优势，但是却会天然使得高质量的合理化学选择空间被压缩，从而拉低整体的输出质量；

以这个分子的反应来进行说明：

compact节省了两倍有余的prompts空间，也适当在合理化学空间当中发生了偏移和侧重；

![q0004_compact_full_case](E:\Python\skills\Rachel\plan\figures\q0004_compact_full_case.png)

![section6_fig2_context_vs_quality](E:\Python\skills\Rachel\plan\figures\section6_fig2_context_vs_quality.png)

也由此，使得我们基于compact设计了渐进式披露的上下文系统来优化化学质量与prompts上下文占用，在面向多步逆合成的过程中，保持在同一个上下文窗口时非常重要的，这可以让llm避开无意义的官能团置换（FGI，例如重复酰卤-酸酐-酯-酰胺-酸衍生的无意义循环）

-------------

### 3     Paroutes120数据benchmark

使用Rachel系统，对于从Paroutes120数据集中提取的120（n1+n5=51+69）条数据进行端到端的全量测试，Paroutes是从PaRoutes的n1\n5各一万条数据库中提取出来的；数据集平均反应步长在3.55；目标复杂度均值大概是 SA=3.012、CS=4.586、BertzCT=1003.58，

目前使用了多种传统方法MCTS+RETRO*与传统模型的对比和Rachel进行效果的对比，但是不得不承认，找到一个统一的评价口径来进行对比Rachel和传统方法是比较困难的，solve rate是只拆分到可买分子，传统方法是使用的社区自建分子中间体库；Rachel是进行的显示可买联网检索审计，solve rate 只能作为一种理论参考，因为很多中间物、terminal尽管可以在反应路线中写出来，但是可能面临着是现制，难以直接检索到成品的，这是合成中经常会发生的一个现象。

<img src="E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\reports\baseline_figures\figures\figure1_solve_rate_n1_only.png" alt="figure1_solve_rate_n1_only"  />![figure1_solve_rate_n5_only](E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\reports\baseline_figures\figures\figure1_solve_rate_n5_only.png)

其中n1\n5数据集更加详细的逆合成深度与从头合成程度详细对比图：![figureC1_common_solved_complexity_compression_n1_only](E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\reports\complexity_analysis\figures\figureC1_common_solved_complexity_compression_n1_only.png)![figureC2_common_solved_route_structure_n1_only](E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\reports\complexity_analysis\figures\figureC2_common_solved_route_structure_n1_only.png)



![figureC1_common_solved_complexity_compression_n5_only](E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\reports\complexity_analysis\figures\figureC1_common_solved_complexity_compression_n5_only.png)![figureC2_common_solved_route_structure_n5_only](E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\reports\complexity_analysis\figures\figureC2_common_solved_route_structure_n5_only.png)

此外，Rachel系统整体给出的反应条件和合成思路也更加详细与完善全面：（都有可视化的文件来进行复盘）![n1_366_groundtruth_vs_rachel_annotated_case_zh](E:\Python\skills\MCTS-RETRO\benchmark\paroutes120\artifacts\reports\case_figures\n1_366_groundtruth_vs_rachel_annotated_case_zh.png)



### 4    使用ReactionT5 forward 模型对于Rachel在paroutes120中的所有反应进行交叉验证

![7433c4c6-a141-4cb4-a85b-cbfed9fa6e1b](E:\Python\skills\Rachel\E7\figures\7433c4c6-a141-4cb4-a85b-cbfed9fa6e1b.png)

使用T5前向预测，对于同一个ps（production），使用T5 forward来交叉预测rs（reaction），只取top10结果来验证Rachel给出的反应是否合理，是否同样在化学合理空间内，结果如上；Rachel给出的结果全部rdkit合理可解析；这一部分是进一步来交叉验证Rachel执行时的逆合成是否合理，能否同样在合理范围内映射到第三方的化学空间中；（top20\30结果同样有，差不多，若是为了数据美观可以都放上去）



### 5     paroutes120的LLM对照效果与实验

使用纯LLM进行同样的分子逆合成分析，在除去smiles不合理、简单可买统计后结果如下，但是未经历化学合理性审计，路线明显存在幻觉和化学错误；但是没有逐一进行单独审计了，llm手写smiles还是挺容易出错的![solve_rate_n1_n5](E:\Python\skills\MCTS-RETRO\LLM\figures\solve_rate_n1_n5.png)

### 6    USPTO190的benchmark

此部分没有ground truth，可买性正在统计中，可以再和其他工作进行横向对比；也已经完成了190各的全部端到端的Rachel从头逆合成分析；但是还是同样的，缺乏统一口径的对比，不同方法的对比存在难度；这部分如果对于具体路线继续优化可以继续提升solverate，但是感觉没有必要了，那就真赛博斗蛐蛐了；所有的实验的端到端实验都是对每一个分子明确使用Rachel进行自动执行（prompts示例：“指明分子smiles目录位置”+“使用Rachel系统对该分子进行化学路线合理、高质量的从头合成”；面向可买的prompts示例：“指明分子smiles目录位置”+“使用Rachel系统对该分子进行化学路线合理、高质量的从头合成”+“完成合成之后再进行可买性审计”），暂时没有采用人工逐轮交互方式（工作量大）

![763ca028-0263-47a9-957d-d3f28681a076](E:\Python\skills\MCTS-RETRO\benchmark\uspto190_online_audit\figures\763ca028-0263-47a9-957d-d3f28681a076.png)

![f128ef7a-c236-47c0-a970-dac9d2f33370](E:\Python\skills\MCTS-RETRO\benchmark\uspto190_online_audit\figures\f128ef7a-c236-47c0-a970-dac9d2f33370.png)



### 7    Rachel的系统组件消融实验

![42707109-12d1-499f-a214-c53199ccd1a2](E:\Python\skills\Rachel\tests\E6\figures\42707109-12d1-499f-a214-c53199ccd1a2.png)

no_validation部分同样暂时没有进行逐一的化学合理性的评价审计，但是根据经验，效果也不会好到那里去.....就是看上去的solve率好看而已。

其中的action包括：断键bond、官能团修饰FGI两种

sandbox则是整个沙盒模块；

validation则是所有的验证门控系统，这部分可以查看相关的harness.md

### 8    可能需要补的其他实验

认为需要对Rachel进行适当的成功、失败案例进行化学深入分析；

额外有的20个ground truth是否需要添加进去进行详细的说明；

## 串联实验的思路

### 1   简短的思路

单步效果（验证形式化空间化学有效与高质量）----单步prompts优化（适配llm基座的上下文管理系统）----paroutes120端到端（验证Rachel系统效果）、paroutes120的纯LLM实验（对比llm效果提升明显）----Rachel ablation（验证组件有效性和必要性，包括action、sandbox、validation）---- ReactionT5交叉验证（交叉验证Rachel的validation效果）----USPTO190benchmark（拓宽横向对比效果与边界）

### 2  采用形式化模板的原因：

1. 相比之下,已有方法比如说模板库,其实存在一个问题,就是对于化学反应的分类不够完善,通常没有办法很好地将同一化学反应家族的反应进行很好的归类与模板提取,会导致模板库的冗余与庞大,进而影响搜索效率；传统模板等库输出也是繁重，不利于后期的上下文管理与优化
2. 已有的深度学习for反应预测的模型,其输出smiles存在很多不可解析的地方,并且环境与配置的依赖比较繁重,设备性能要求有门槛,速度也是问题

### 3. 上下文中compact与full的压缩

1. 这一部分其实可以详细见harness中的渐进式披露，rdkit工具提取出来的信息是需要进行自然语言关键提取与转换才能够确保不丢失化学质量的前体下压缩上下文空间的
2. 其实更多的是Rachel构建过程中逐步试探出来的解决办法

### 4. solve rate的价值

1. 只能说指标很理想，但是统计口径很难统一，既有文章叙事的需要，也有现实因素的影响，这一部分仍然建议作为参考；
2. 化学合理性以及完整的路线是相关的工作可能比较没有那么兼顾与侧重的地方

### 5. Rachel里面LLM如何进行分子编辑与断键（详见harness）

1. 首先基于形式化模板+渐进式上下文获取列表
2. LLM作为决策判断与选取，若不满意，其自行进入沙盒，提出方案，若是断键方案，先走预设的断键思路，不满意再进行受限图编辑；若是FGI方案，也是先走预设FGI置换，不满意再进行受限图编辑；
3. 所有沙盒选择都会经过validation

### 6. Rachel如何验证

其实也可以详细见Harness，从原子数、骨架、官能团变化、断键前后重原子数、FGI冲突性等多方面进行了约束和限制。

