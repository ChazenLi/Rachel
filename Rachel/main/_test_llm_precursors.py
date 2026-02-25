#!/usr/bin/env python
"""
LLM 自主提前体 — 深度验证
==========================
模拟 LLM 作为化学家，在以下场景中自主设计前体:

场景 1: 酰胺键拆分 — 模板给酰胺偶联(羧酸+胺)，LLM 改用酰氯+胺（更高效）
场景 2: 卤素优化 — 模板给溴代物，LLM 改用碘代物（Pd偶联更活泼）
场景 3: 保护基策略 — LLM 主动给胺加 Boc 保护再做其他反应
场景 4: 模板失灵 — 复杂分子无模板匹配，LLM 自主设计逆合成路线
场景 5: 多轮试错 — LLM 提了不好的前体，验证失败，修正后重试
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import RDLogger, Chem
RDLogger.DisableLog("rdApp.*")

from Rachel.main import (
    RetrosynthesisOrchestrator,
    SandboxResult,
    LLMDecision,
    generate_forward_report,
    get_terminal_list,
)


def pp(label, sandbox: SandboxResult):
    """Pretty-print sandbox result."""
    sd = sandbox.to_dict()
    print(f"\n  [{label}]")
    print(f"    success: {sd['success']}")
    if sd.get("error"):
        print(f"    error: {sd['error']}")
        return
    if sd.get("precursors"):
        print(f"    precursors: {' + '.join(sd['precursors'])}")
    for pd in sd.get("precursor_details", []):
        tag = "✓term" if pd["is_terminal"] else "…pend"
        seen = " ⚠SEEN" if pd.get("already_seen") else ""
        print(f"      {tag} {pd['smiles'][:55]}  "
              f"CS={pd['cs_score']:.2f} heavy={pd['heavy_atoms']}{seen}")
    fv = sd.get("forward_validation", {})
    if fv:
        icon = "✓" if fv.get("pass") else "✗"
        print(f"    forward: {icon} feasibility={fv.get('feasibility_score', '?')}")
        if fv.get("hard_fail_reasons"):
            print(f"    hard_fail: {fv['hard_fail_reasons']}")
    ab = sd.get("atom_balance", {})
    if ab:
        print(f"    atom_balance: {'✓' if ab.get('balanced') else '✗'} {ab.get('message', '')}")
    if sd.get("cycle_warnings"):
        print(f"    ⚠ cycle: {sd['cycle_warnings']}")


# ═══════════════════════════════════════════════════════════════
# 场景 1: 酰胺键拆分 — 酰氯 vs 羧酸
# ═══════════════════════════════════════════════════════════════

def scenario_1_amide_acyl_chloride():
    """
    目标: N-苯基苯甲酰胺  PhCONHPh
    模板会给: PhCOOH + PhNH2 (酰胺偶联)
    LLM 化学直觉: 酰氯 PhCOCl + PhNH2 更直接，不需要偶联试剂
    """
    print("=" * 70)
    print("场景 1: 酰胺键 — 酰氯 vs 羧酸 vs 酸酐")
    print("=" * 70)

    target = "O=C(Nc1ccccc1)c1ccccc1"  # N-苯基苯甲酰胺
    orch = RetrosynthesisOrchestrator(
        target, target_name="N-Phenylbenzamide",
        terminal_cs_threshold=1.2,
    )

    ctx = orch.prepare_next()
    print(f"目标: {ctx.smiles}")
    print(f"CS={ctx.cs_score:.2f}  tier={ctx.decision_tier}")

    # 先看模板给了什么
    compact = ctx.to_dict()
    print(f"\n模板提供 {compact.get('n_bonds', 0)} 个键位:")
    for bs in compact.get("bond_summary", []):
        print(f"  [{bs['bond_idx']}] atoms={bs['atoms']} "
              f"score={bs['heuristic_score']:.3f} "
              f"types={bs['reaction_types']}")

    # 看看模板的酰胺偶联方案
    print("\n--- 模板方案 (explore_bond) ---")
    for i in range(min(3, compact.get("n_bonds", 0))):
        detail = orch.explore_bond(i)
        for j, alt in enumerate(detail.get("alternatives", [])[:2]):
            print(f"  bond[{i}] alt[{j}]: {alt['template']}")
            print(f"    → {' + '.join(alt['precursors'])}")

    # === LLM 思考 ===
    # "模板给的是 Amide Coupling (羧酸+胺+偶联试剂)，
    #  但我知道酰氯+胺是更经典的 Schotten-Baumann 反应，
    #  不需要 EDC/HATU，条件更温和。让我试试。"

    print("\n--- LLM 自主方案 A: 酰氯 + 胺 (Schotten-Baumann) ---")
    sb_a = orch.try_precursors(
        ["O=C(Cl)c1ccccc1", "Nc1ccccc1"],
        reaction_type="Schotten-Baumann acylation",
    )
    pp("酰氯+苯胺", sb_a)

    # === LLM 继续思考 ===
    # "酰氯方案验证通过了。但让我也试试酸酐方案，
    #  苯甲酸酐 + 苯胺也是经典方法。"

    print("\n--- LLM 自主方案 B: 酸酐 + 胺 ---")
    sb_b = orch.try_precursors(
        ["O=C(OC(=O)c1ccccc1)c1ccccc1", "Nc1ccccc1"],
        reaction_type="anhydride aminolysis",
    )
    pp("苯甲酸酐+苯胺", sb_b)

    # === LLM 决策 ===
    # "两个方案都能过验证。酰氯方案前体更简单(CS更低)，
    #  而且酰氯+胺是本科有机化学最基础的反应，选它。"

    best = sb_a if sb_a.success else sb_b
    if best.success:
        result = orch.commit_decision(
            precursor_override=best.precursors,
            reaction_type=best.reaction_type,
            llm_decision=LLMDecision(
                selection_reasoning="酰氯+胺 Schotten-Baumann 反应，"
                    "不需偶联试剂，条件温和，前体简单易得",
                confidence="high",
                rejected_alternatives=[
                    {"method": "amide coupling (HATU)", "reason": "需要昂贵偶联试剂"},
                    {"method": "anhydride aminolysis", "reason": "酸酐前体更复杂"},
                ],
            ),
        )
        print(f"\n  ✓ commit: {result.success}")
        print(f"    terminal: {result.new_terminal}")
        print(f"    pending: {result.new_pending}")

    # 清理
    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None: break
        orch.accept_terminal(reason="cleanup")

    orch.finalize("LLM chose Schotten-Baumann over template amide coupling")
    print(f"\n{generate_forward_report(orch.get_tree())}")
    print("\n✓ 场景 1 完成")


# ═══════════════════════════════════════════════════════════════
# 场景 2: 卤素优化 — 溴→碘 (Pd偶联更活泼)
# ═══════════════════════════════════════════════════════════════

def scenario_2_halogen_swap():
    """
    目标: 4-甲基联苯  4-Me-C6H4-C6H5
    模板可能给: ArBr + PhB(OH)2 (Suzuki)
    LLM 化学直觉: ArI 在 Pd 催化偶联中反应活性更高，
    尤其对于缺电子底物，碘代物转化率更好。
    """
    print("\n" + "=" * 70)
    print("场景 2: 卤素优化 — Suzuki 偶联中溴→碘")
    print("=" * 70)

    target = "Cc1ccc(-c2ccccc2)cc1"  # 4-甲基联苯
    orch = RetrosynthesisOrchestrator(
        target, target_name="4-Methylbiphenyl",
        terminal_cs_threshold=1.2,
    )

    ctx = orch.prepare_next()
    print(f"目标: {ctx.smiles}  CS={ctx.cs_score:.2f}")

    # 看模板方案
    compact = ctx.to_dict()
    print(f"\n模板键位概览:")
    for bs in compact.get("bond_summary", []):
        print(f"  [{bs['bond_idx']}] atoms={bs['atoms']} "
              f"score={bs['heuristic_score']:.3f} types={bs['reaction_types']}")

    # 模板的 Suzuki 方案（通常给 ArBr）
    print("\n--- 模板 Suzuki 方案 ---")
    detail = orch.explore_bond(0)
    for j, alt in enumerate(detail.get("alternatives", [])[:3]):
        print(f"  alt[{j}]: {alt['template']}")
        print(f"    → {' + '.join(alt['precursors'])}")

    # 模板沙盒试断
    sb_template = orch.try_disconnection(0, 0)
    pp("模板方案", sb_template)

    # === LLM 思考 ===
    # "模板给了 ArBr + PhB(OH)2。但我知道 ArI 在 Suzuki 偶联中
    #  氧化加成更快，尤其室温条件下碘代物就能反应。
    #  让我把溴换成碘试试。"

    print("\n--- LLM 方案: ArI + PhB(OH)2 (碘代物更活泼) ---")
    sb_iodo = orch.try_precursors(
        ["Ic1ccc(C)cc1", "OB(O)c1ccccc1"],
        reaction_type="Suzuki coupling (ArI)",
    )
    pp("碘代物Suzuki", sb_iodo)

    # === LLM 也试试反向极性 ===
    print("\n--- LLM 方案: ArB(OH)2 + PhBr (反向极性) ---")
    sb_reverse = orch.try_precursors(
        ["OB(O)c1ccc(C)cc1", "Brc1ccccc1"],
        reaction_type="Suzuki coupling (reverse)",
    )
    pp("反向Suzuki", sb_reverse)

    # 选最佳方案 commit
    best = sb_iodo if sb_iodo.success else sb_template
    if best.success:
        result = orch.commit_decision(
            precursor_override=best.precursors,
            reaction_type=best.reaction_type,
            llm_decision=LLMDecision(
                selection_reasoning="ArI 在 Pd 催化 Suzuki 偶联中氧化加成更快，"
                    "室温即可反应，转化率更高",
                confidence="high",
            ),
        )
        print(f"\n  ✓ commit: {result.success}")

    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None: break
        orch.accept_terminal(reason="cleanup")

    orch.finalize("LLM optimized halogen: Br→I for Suzuki")
    print(f"\n{generate_forward_report(orch.get_tree())}")
    print("\n✓ 场景 2 完成")


# ═══════════════════════════════════════════════════════════════
# 场景 3: 保护基策略 — LLM 主动加 Boc 保护
# ═══════════════════════════════════════════════════════════════

def scenario_3_protection_strategy():
    """
    目标: 4-氨基苯甲酸甲酯  H2N-C6H4-COOMe
    问题: 如果直接做酯化，游离胺基会干扰
    LLM 策略: 先拆酯键，但给胺加 Boc 保护
    """
    print("\n" + "=" * 70)
    print("场景 3: 保护基策略 — Boc 保护胺基")
    print("=" * 70)

    target = "COC(=O)c1ccc(N)cc1"  # 4-氨基苯甲酸甲酯
    orch = RetrosynthesisOrchestrator(
        target, target_name="Methyl 4-aminobenzoate",
        terminal_cs_threshold=1.2,
    )

    ctx = orch.prepare_next()
    print(f"目标: {ctx.smiles}  CS={ctx.cs_score:.2f}")

    compact = ctx.to_dict()
    print(f"\n模板键位: {compact.get('n_bonds', 0)} 个")
    for bs in compact.get("bond_summary", []):
        print(f"  [{bs['bond_idx']}] atoms={bs['atoms']} "
              f"score={bs['heuristic_score']:.3f} types={bs['reaction_types']}")

    # === LLM 思考 ===
    # "这个分子有游离胺基和酯基。模板可能直接拆酯键给
    #  4-氨基苯甲酸 + MeOH。但在实际合成中，游离胺基
    #  在很多反应条件下会出问题。
    #
    #  更好的策略: 先认为这个分子来自 Boc 脱保护，
    #  即 Boc-NH-C6H4-COOMe → H2N-C6H4-COOMe
    #  这样上游合成时胺基是被保护的。"

    print("\n--- LLM 方案 A: 直接拆酯 (模板思路) ---")
    # 先看模板怎么拆
    sb_ester = orch.try_precursors(
        ["OC(=O)c1ccc(N)cc1", "CO"],
        reaction_type="Fischer esterification",
    )
    pp("直接拆酯", sb_ester)

    print("\n--- LLM 方案 B: Boc 脱保护 (保护基策略) ---")
    sb_boc = orch.try_precursors(
        ["COC(=O)c1ccc(NC(=O)OC(C)(C)C)cc1"],
        reaction_type="Boc deprotection (TFA)",
    )
    pp("Boc脱保护", sb_boc)

    print("\n--- LLM 方案 C: 硝基还原 (从硝基苯甲酸甲酯) ---")
    sb_nitro = orch.try_precursors(
        ["COC(=O)c1ccc([N+](=O)[O-])cc1"],
        reaction_type="nitro reduction (H2/Pd)",
    )
    pp("硝基还原", sb_nitro)

    # === LLM 决策 ===
    # "三个方案都验证通过了:
    #  A) 直接拆酯 — 简单但游离胺基在酯化条件下可能有副反应
    #  B) Boc脱保护 — 保护基策略，上游合成更安全
    #  C) 硝基还原 — 经典方法，硝基苯甲酸甲酯是商品化试剂
    #  选 C，因为硝基苯甲酸甲酯便宜易得，还原条件温和。"

    best = sb_nitro if sb_nitro.success else sb_boc
    if best.success:
        result = orch.commit_decision(
            precursor_override=best.precursors,
            reaction_type=best.reaction_type,
            llm_decision=LLMDecision(
                selection_reasoning="硝基还原(H2/Pd-C)条件温和、选择性好，"
                    "硝基苯甲酸甲酯是廉价商品化试剂",
                confidence="high",
                rejected_alternatives=[
                    {"method": "direct esterification",
                     "reason": "游离胺基在酸性酯化条件下会被质子化"},
                    {"method": "Boc deprotection",
                     "reason": "多一步保护/脱保护，不如硝基还原直接"},
                ],
                protection_needed=False,
            ),
        )
        print(f"\n  ✓ commit: {result.success}")
        print(f"    terminal: {result.new_terminal}")
        print(f"    pending: {result.new_pending}")

    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None: break
        orch.accept_terminal(reason="cleanup")

    orch.finalize("LLM chose nitro reduction over direct esterification")
    print(f"\n{generate_forward_report(orch.get_tree())}")
    print("\n✓ 场景 3 完成")


# ═══════════════════════════════════════════════════════════════
# 场景 4: 模板失灵 — LLM 自主设计路线
# ═══════════════════════════════════════════════════════════════

def scenario_4_template_failure():
    """
    目标: 一个模板可能处理不好的分子
    2-苯基吲哚  — Fischer 吲哚合成是经典方法
    LLM 知道: 苯肼 + 苯乙酮 → 2-苯基吲哚 (Fischer indole)
    """
    print("\n" + "=" * 70)
    print("场景 4: 模板失灵 — LLM 自主设计 Fischer 吲哚合成")
    print("=" * 70)

    target = "c1ccc(-c2cc3ccccc3[nH]2)cc1"  # 2-苯基吲哚
    orch = RetrosynthesisOrchestrator(
        target, target_name="2-Phenylindole",
        terminal_cs_threshold=1.2,
    )

    ctx = orch.prepare_next()
    print(f"目标: {ctx.smiles}  CS={ctx.cs_score:.2f}")

    compact = ctx.to_dict()
    n_bonds = compact.get("n_bonds", 0)
    print(f"模板键位: {n_bonds} 个")

    # 看看模板能给什么
    if n_bonds > 0:
        for i in range(min(3, n_bonds)):
            detail = orch.explore_bond(i)
            alts = detail.get("alternatives", [])
            if alts:
                print(f"  bond[{i}]: {alts[0]['template']}")
                print(f"    → {' + '.join(alts[0]['precursors'])}")
            else:
                print(f"  bond[{i}]: atoms={detail['atoms']} — 无可用模板")

    # === LLM 思考 ===
    # "2-苯基吲哚的经典合成是 Fischer 吲哚合成:
    #  苯肼 + 苯乙酮 → 在酸催化下经 [3,3]-σ 迁移得到吲哚。
    #  这是有机化学教科书级别的反应，模板可能没有覆盖。
    #  让我直接提出这个方案。"

    print("\n--- LLM 方案: Fischer 吲哚合成 ---")
    print("  苯肼 + 苯乙酮 → 2-苯基吲哚 (酸催化, [3,3]-σ迁移)")
    sb_fischer = orch.try_precursors(
        ["NNc1ccccc1", "CC(=O)c1ccccc1"],
        reaction_type="Fischer indole synthesis",
    )
    pp("Fischer吲哚", sb_fischer)

    # === LLM 也试试 Larock 吲哚合成 ===
    print("\n--- LLM 方案: Larock 吲哚合成 ---")
    print("  2-碘苯胺 + 苯乙炔 → 2-苯基吲哚 (Pd催化)")
    sb_larock = orch.try_precursors(
        ["Nc1ccccc1I", "C#Cc1ccccc1"],
        reaction_type="Larock indole synthesis",
    )
    pp("Larock吲哚", sb_larock)

    # === LLM 决策 ===
    # "Fischer 方案前体更简单(苯肼+苯乙酮都是常见试剂)，
    #  不需要 Pd 催化剂。选 Fischer。"

    best = sb_fischer if sb_fischer.success else sb_larock
    if best.success:
        result = orch.commit_decision(
            precursor_override=best.precursors,
            reaction_type=best.reaction_type,
            llm_decision=LLMDecision(
                selection_reasoning="Fischer 吲哚合成是经典 name reaction，"
                    "苯肼+苯乙酮均为廉价商品化试剂，酸催化条件简单",
                confidence="high",
                rejected_alternatives=[
                    {"method": "Larock indole", "reason": "需要 Pd 催化剂和碘代物"},
                ],
            ),
        )
        print(f"\n  ✓ commit: {result.success}")
        print(f"    terminal: {result.new_terminal}")
        print(f"    pending: {result.new_pending}")
    else:
        print("\n  两个方案都失败了，skip")
        orch.skip_current("Fischer and Larock both failed")

    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None: break
        orch.accept_terminal(reason="cleanup")

    orch.finalize("LLM designed Fischer indole synthesis")
    print(f"\n{generate_forward_report(orch.get_tree())}")
    print("\n✓ 场景 4 完成")


# ═══════════════════════════════════════════════════════════════
# 场景 5: 多轮试错 — LLM 反复修正前体
# ═══════════════════════════════════════════════════════════════

def scenario_5_iterative_refinement():
    """
    目标: 4-氨基-3-硝基苯甲酸  (有胺基+硝基+羧酸，官能团兼容性复杂)
    LLM 需要多轮试错来找到合适的前体组合。
    """
    print("\n" + "=" * 70)
    print("场景 5: 多轮试错 — 反复修正前体")
    print("=" * 70)

    target = "Nc1ccc(C(=O)O)cc1[N+](=O)[O-]"  # 4-氨基-3-硝基苯甲酸
    orch = RetrosynthesisOrchestrator(
        target, target_name="4-Amino-3-nitrobenzoic acid",
        terminal_cs_threshold=1.3,
    )

    ctx = orch.prepare_next()
    print(f"目标: {ctx.smiles}  CS={ctx.cs_score:.2f}")

    compact = ctx.to_dict()
    print(f"模板键位: {compact.get('n_bonds', 0)} 个")
    print(f"FGI选项: {compact.get('n_fgi', 0)} 个")

    # === LLM 第一轮: 天真的尝试 ===
    print("\n--- 第 1 轮: 直接拆硝基 (硝化反应逆合成) ---")
    print("  思路: 4-氨基苯甲酸 + 硝化 → 产物")
    sb1 = orch.try_precursors(
        ["Nc1ccc(C(=O)O)cc1"],  # 4-氨基苯甲酸
        reaction_type="nitration",
    )
    pp("硝化(单前体)", sb1)

    # === LLM 第二轮: 换个思路 ===
    print("\n--- 第 2 轮: 先还原硝基 → 从二硝基出发 ---")
    print("  思路: 3,4-二硝基苯甲酸 选择性还原一个硝基")
    sb2 = orch.try_precursors(
        ["O=C(O)c1ccc([N+](=O)[O-])c([N+](=O)[O-])c1"],
        reaction_type="selective nitro reduction",
    )
    pp("选择性还原", sb2)

    # === LLM 第三轮: 经典方法 ===
    print("\n--- 第 3 轮: 从 4-硝基苯甲酸出发 ---")
    print("  思路: 4-硝基苯甲酸 → 硝化得到 3,4-二硝基 → 选择性还原")
    print("  但这一步只做: 二硝基苯甲酸 → 选择性还原 → 产物")
    print("  不对，让我重新想...")

    print("\n--- 第 3 轮(修正): 从 4-氯-3-硝基苯甲酸 + 氨 ---")
    print("  思路: SNAr 用氨替换氯 (硝基邻位活化)")
    sb3 = orch.try_precursors(
        ["O=C(O)c1ccc(Cl)c([N+](=O)[O-])c1", "N"],
        reaction_type="SNAr amination",
    )
    pp("SNAr氨化", sb3)

    # === LLM 第四轮: 最终方案 ===
    print("\n--- 第 4 轮: 从 4-氨基苯甲酸 硝化 (重新审视) ---")
    print("  思路: 虽然第1轮试过，但硝化是单前体反应，")
    print("  让我换个角度 — 从 3-硝基-4-氯苯甲酸 + NH3 更合理")

    # 比较所有成功方案
    print("\n=== LLM 决策总结 ===")
    results = [
        ("硝化(单前体)", sb1),
        ("选择性还原", sb2),
        ("SNAr氨化", sb3),
    ]
    for name, sb in results:
        if sb.success:
            fv = sb.forward_validation or {}
            score = fv.get("assessment", {}).get("feasibility_score", "N/A")
            prec = " + ".join(sb.precursors)
            print(f"  {name}: feasibility={score}  → {prec[:60]}")
        else:
            print(f"  {name}: FAILED ({sb.error})")

    # 选最佳
    successful = [(n, s) for n, s in results if s.success]
    if successful:
        best_name, best_sb = successful[0]
        # 如果有多个成功的，选 feasibility 最高的
        for n, s in successful[1:]:
            fv_best = (best_sb.forward_validation or {}).get(
                "assessment", {}).get("feasibility_score", 0)
            fv_curr = (s.forward_validation or {}).get(
                "assessment", {}).get("feasibility_score", 0)
            if fv_curr and fv_curr > (fv_best or 0):
                best_name, best_sb = n, s

        print(f"\n  → 选择: {best_name}")
        result = orch.commit_decision(
            precursor_override=best_sb.precursors,
            reaction_type=best_sb.reaction_type,
            llm_decision=LLMDecision(
                selection_reasoning=f"经过 {len(results)} 轮试错，"
                    f"选择 {best_name} 方案",
                confidence="medium",
            ),
        )
        print(f"  commit: {result.success}")
    else:
        print("\n  所有方案都失败了")
        orch.skip_current("all attempts failed")

    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None: break
        orch.accept_terminal(reason="cleanup")

    orch.finalize("LLM iterative refinement")
    print(f"\n{generate_forward_report(orch.get_tree())}")
    print("\n✓ 场景 5 完成")


# ═══════════════════════════════════════════════════════════════
# 场景 6: 多步自主路线 — LLM 全程自主设计
# ═══════════════════════════════════════════════════════════════

def scenario_6_full_autonomous():
    """
    目标: 对乙酰氨基酚 (扑热息痛)  4-HO-C6H4-NHCOCH3
    LLM 全程自主设计路线，不依赖模板。
    经典工业路线: 对氨基酚 + 乙酸酐 → 扑热息痛
    """
    print("\n" + "=" * 70)
    print("场景 6: 全程自主 — 扑热息痛合成路线设计")
    print("=" * 70)

    target = "CC(=O)Nc1ccc(O)cc1"  # 对乙酰氨基酚
    orch = RetrosynthesisOrchestrator(
        target, target_name="Paracetamol",
        terminal_cs_threshold=1.2,
    )

    # === Step 1: 拆酰胺键 ===
    ctx = orch.prepare_next()
    print(f"\n[Step 1] 目标: {ctx.smiles}  CS={ctx.cs_score:.2f}")

    # LLM: "扑热息痛的酰胺键是最明显的断键位。
    #  工业上用对氨基酚 + 乙酸酐，我来验证。"

    print("  LLM思考: 拆酰胺键 — 对氨基酚 + 乙酸酐")
    sb1 = orch.try_precursors(
        ["Nc1ccc(O)cc1", "CC(=O)OC(C)=O"],
        reaction_type="acetylation (anhydride)",
    )
    pp("乙酸酐法", sb1)

    # 也试试乙酰氯
    print("  LLM思考: 或者用乙酰氯?")
    sb1b = orch.try_precursors(
        ["Nc1ccc(O)cc1", "CC(=O)Cl"],
        reaction_type="acetylation (acyl chloride)",
    )
    pp("乙酰氯法", sb1b)

    # 选乙酸酐（工业标准）
    best1 = sb1 if sb1.success else sb1b
    if best1.success:
        result = orch.commit_decision(
            precursor_override=best1.precursors,
            reaction_type=best1.reaction_type,
            llm_decision=LLMDecision(
                selection_reasoning="乙酸酐法是工业标准路线，"
                    "副产物乙酸可回收，原子经济性好",
                confidence="high",
            ),
        )
        print(f"\n  ✓ commit step 1: {result.success}")
        for t in result.new_terminal:
            print(f"    ✓ terminal: {t}")
        for p in result.new_pending:
            print(f"    … pending: {p}")

    # === Step 2: 对氨基酚怎么来? ===
    ctx2 = orch.prepare_next()
    if ctx2 and ctx2.decision_tier != "quick_pass":
        print(f"\n[Step 2] 当前: {ctx2.smiles}  CS={ctx2.cs_score:.2f}")

        # LLM: "对氨基酚的工业合成:
        #  方法A: 对硝基酚 → 催化加氢还原
        #  方法B: 苯酚 → 硝化 → 还原 (但选择性问题)
        #  选方法A，对硝基酚是商品化试剂。"

        print("  LLM思考: 对硝基酚 催化加氢还原")
        sb2 = orch.try_precursors(
            ["O=[N+]([O-])c1ccc(O)cc1"],
            reaction_type="catalytic hydrogenation",
        )
        pp("硝基还原", sb2)

        if sb2.success:
            result2 = orch.commit_decision(
                precursor_override=sb2.precursors,
                reaction_type=sb2.reaction_type,
                llm_decision=LLMDecision(
                    selection_reasoning="对硝基酚是廉价商品化试剂，"
                        "H2/Pd-C 还原条件温和、选择性好",
                    confidence="high",
                ),
            )
            print(f"\n  ✓ commit step 2: {result2.success}")
            for t in result2.new_terminal:
                print(f"    ✓ terminal: {t}")
            for p in result2.new_pending:
                print(f"    … pending: {p}")
    elif ctx2:
        print(f"\n[Step 2] {ctx2.smiles} → quick_pass terminal")
        orch.accept_terminal(reason="quick_pass")

    # 清理剩余
    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None: break
        if c.decision_tier == "quick_pass":
            orch.accept_terminal(reason="quick_pass")
        else:
            orch.skip_current("cleanup")

    orch.finalize("LLM autonomous: industrial paracetamol route")

    tree = orch.get_tree()
    print(f"\n{generate_forward_report(tree)}")
    print(f"\n{tree.print_tree()}")

    terminals = get_terminal_list(tree)
    print(f"\n起始原料 ({len(terminals)} 种):")
    for t in terminals:
        print(f"  {t['smiles']}  CS={t.get('cs_score', 0):.2f}")

    print("\n✓ 场景 6 完成")


# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    scenario_1_amide_acyl_chloride()
    scenario_2_halogen_swap()
    scenario_3_protection_strategy()
    scenario_4_template_failure()
    scenario_5_iterative_refinement()
    scenario_6_full_autonomous()

    print("\n" + "=" * 70)
    print("✓ 全部 6 个场景测试完成")
    print("=" * 70)
