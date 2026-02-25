#!/usr/bin/env python
"""测试沙盒试探 + 分层上下文 — 模拟 LLM 交互流程。

演示完整的 LLM 决策循环:
  prepare_next()       → 精简摘要
  explore_bond(idx)    → 按需展开
  try_disconnection()  → 沙盒试断
  try_precursors()     → LLM 自提前体验证
  commit_decision()    → 正式写入
"""
import sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.main import (
    RetrosynthesisOrchestrator,
    LLMDecision,
    generate_forward_report,
)


def test_layered_context():
    """测试 1: 分层上下文 — compact vs full。"""
    print("=" * 60)
    print("测试 1: 分层上下文")
    print("=" * 60)

    orch = RetrosynthesisOrchestrator(
        "CC(=O)Oc1ccccc1C(=O)O",
        target_name="Aspirin",
        terminal_cs_threshold=1.2,
    )

    ctx = orch.prepare_next()
    assert ctx is not None
    assert ctx.decision_tier == "standard"

    # compact 模式（默认）
    compact = ctx.to_dict()
    print(f"\n[compact] keys: {sorted(compact.keys())}")
    print(f"  bond_summary: {len(compact.get('bond_summary', []))} 个键位概览")
    assert "bond_summary" in compact
    assert "disconnectable_bonds" not in compact  # compact 不含完整方案
    for bs in compact.get("bond_summary", [])[:2]:
        print(f"    [{bs['bond_idx']}] atoms={bs['atoms']} "
              f"score={bs['heuristic_score']:.3f} "
              f"n_alt={bs['n_alternatives']} "
              f"types={bs['reaction_types']}")

    # full 模式
    full = ctx.to_dict(detail="full")
    print(f"\n[full] keys: {sorted(full.keys())}")
    assert "disconnectable_bonds" in full
    n_alts = sum(len(b.get("alternatives", []))
                 for b in full.get("disconnectable_bonds", []))
    print(f"  disconnectable_bonds: {len(full['disconnectable_bonds'])} 个键位, "
          f"共 {n_alts} 个方案")

    # compact 文本
    text = ctx.to_text()
    print(f"\n[compact text]:\n{text}")

    # 清理
    orch.skip_current("test done")
    print("\n✓ 分层上下文测试通过")


def test_explore_and_sandbox():
    """测试 2: explore + sandbox 试探。"""
    print("\n" + "=" * 60)
    print("测试 2: explore + sandbox 试探")
    print("=" * 60)

    orch = RetrosynthesisOrchestrator(
        "CC(=O)Oc1ccccc1C(=O)O",
        target_name="Aspirin",
        terminal_cs_threshold=1.2,
    )

    ctx = orch.prepare_next()
    assert ctx is not None

    # ── explore_bond ──
    print("\n--- explore_bond(0) ---")
    detail = orch.explore_bond(0)
    print(f"  atoms={detail['atoms']}  score={detail['heuristic_score']:.3f}")
    alts = detail.get("alternatives", [])
    print(f"  alternatives: {len(alts)}")
    for i, alt in enumerate(alts[:3]):
        print(f"    [{i}] {alt['template']}")
        print(f"        → {' + '.join(alt['precursors'])}")

    # ── try_disconnection ──
    print("\n--- try_disconnection(0, 0) ---")
    sandbox = orch.try_disconnection(0, 0)
    sd = sandbox.to_dict()
    print(f"  success={sd['success']}")
    if sd["success"]:
        print(f"  precursors: {sd['precursors']}")
        for pd in sd.get("precursor_details", []):
            print(f"    {pd['smiles'][:50]}  CS={pd['cs_score']:.2f}  "
                  f"terminal={pd['is_terminal']}  heavy={pd['heavy_atoms']}")
        fv = sd.get("forward_validation", {})
        print(f"  forward: pass={fv.get('pass')}  "
              f"feasibility={fv.get('feasibility_score')}")
        ab = sd.get("atom_balance", {})
        print(f"  atom_balance: {ab.get('balanced')}")

    # ── 试另一个键位 ──
    print("\n--- try_disconnection(1, 0) ---")
    sandbox2 = orch.try_disconnection(1, 0)
    sd2 = sandbox2.to_dict()
    print(f"  success={sd2['success']}")
    if sd2["success"]:
        print(f"  precursors: {sd2['precursors']}")

    # ── try_precursors: LLM 自提前体 ──
    print("\n--- try_precursors (LLM 自提) ---")
    sandbox3 = orch.try_precursors(
        ["c1ccccc1O", "CC(=O)Cl"],
        reaction_type="acylation",
    )
    sd3 = sandbox3.to_dict()
    print(f"  success={sd3['success']}")
    if sd3["success"]:
        print(f"  precursors: {sd3['precursors']}")
        for pd in sd3.get("precursor_details", []):
            print(f"    {pd['smiles']}  CS={pd['cs_score']:.2f}  "
                  f"terminal={pd['is_terminal']}")
        fv = sd3.get("forward_validation", {})
        print(f"  forward: pass={fv.get('pass')}  "
              f"feasibility={fv.get('feasibility_score')}")

    # ── 无效 SMILES 测试 ──
    print("\n--- try_precursors (invalid SMILES) ---")
    sandbox4 = orch.try_precursors(["INVALID_SMILES", "CC"])
    print(f"  success={sandbox4.success}  error={sandbox4.error}")

    # ── 最终 commit ──
    print("\n--- commit_decision (选择方案 0,0) ---")
    if alts:
        alt = alts[0]
        result = orch.commit_decision(
            bond=tuple(detail["atoms"]),
            reaction_type=alt["template"].split("(")[0].strip(),
            template_id=alt.get("template_id", ""),
            template_name=alt["template"],
            llm_decision=LLMDecision(
                selection_reasoning="沙盒试探后选择 heuristic 最高方案",
                confidence="high",
            ),
        )
        print(f"  success={result.success}")
        if result.success:
            print(f"  new_terminal: {result.new_terminal}")
            print(f"  new_pending: {result.new_pending}")

    # 处理剩余
    while not orch.is_complete():
        c = orch.prepare_next()
        if c is None:
            break
        orch.accept_terminal(reason="test cleanup")

    print(f"\n  Final: {orch.get_status()}")
    print("\n✓ explore + sandbox 测试通过")


def test_explore_fgi():
    """测试 3: FGI 探索 + 沙盒。"""
    print("\n" + "=" * 60)
    print("测试 3: FGI explore + sandbox")
    print("=" * 60)

    # 用一个有 FGI 选项的分子
    orch = RetrosynthesisOrchestrator(
        "OB(O)c1ccccc1",
        target_name="Phenylboronic acid",
        terminal_cs_threshold=1.0,
    )

    ctx = orch.prepare_next()
    if ctx is None or ctx.decision_tier == "quick_pass":
        print("  分子太简单，跳过 FGI 测试")
        orch.accept_terminal(reason="too simple")
    else:
        fgi_info = orch.explore_fgi()
        n_fgi = fgi_info.get("n_fgi", 0)
        print(f"  FGI options: {n_fgi}")
        if n_fgi > 0:
            for i, f in enumerate(fgi_info["fgi_options"][:3]):
                print(f"    [{i}] {f['template']} → {' + '.join(f['precursors'])}")

            sandbox = orch.try_fgi(0)
            sd = sandbox.to_dict()
            print(f"\n  try_fgi(0): success={sd['success']}")
            if sd["success"]:
                print(f"    precursors: {sd['precursors']}")
        else:
            print("  无 FGI 选项")
        orch.skip_current("test done")

    print("\n✓ FGI 测试通过")


def test_full_imatinib_with_sandbox():
    """测试 4: Imatinib 完整流程 — 模拟 LLM 用沙盒试探后决策。"""
    print("\n" + "=" * 60)
    print("测试 4: Imatinib — 沙盒辅助决策")
    print("=" * 60)

    orch = RetrosynthesisOrchestrator(
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        target_name="Imatinib",
        terminal_cs_threshold=1.5,
        max_depth=5,
        max_steps=20,
    )

    step = 0
    while not orch.is_complete():
        ctx = orch.prepare_next()
        if ctx is None:
            break

        if ctx.decision_tier == "quick_pass":
            print(f"  ✓ terminal: {ctx.smiles[:50]}  CS={ctx.cs_score:.2f}")
            orch.accept_terminal(reason="quick_pass")
            continue

        step += 1
        print(f"\n[Step {step}] {ctx.smiles[:60]}  CS={ctx.cs_score:.2f}")

        # 模拟 LLM: 先看 compact 摘要
        compact = ctx.to_dict()
        bonds = compact.get("bond_summary", [])
        print(f"  compact: {len(bonds)} 键位概览")

        if not bonds:
            orch.skip_current("no bonds")
            continue

        # 模拟 LLM: 对 top-1 键位展开详情
        detail = orch.explore_bond(0)
        alts = detail.get("alternatives", [])
        print(f"  explore_bond(0): {len(alts)} 方案")

        if not alts:
            orch.skip_current("no alternatives")
            continue

        # 模拟 LLM: 沙盒试断 top-1
        sandbox = orch.try_disconnection(0, 0)
        if sandbox.success:
            prec_str = " + ".join(sandbox.precursors)
            print(f"  sandbox: {sandbox.template_name}")
            print(f"    → {prec_str[:70]}")
            fv = sandbox.forward_validation
            if fv:
                score = fv.get("assessment", {}).get("feasibility_score", 0)
                print(f"    feasibility={score:.3f}")

            # 满意，commit
            result = orch.commit_decision(
                bond=tuple(detail["atoms"]),
                reaction_type=alts[0].get("template", "").split("(")[0].strip(),
                template_id=alts[0].get("template_id", ""),
                template_name=alts[0].get("template", ""),
                llm_decision=LLMDecision(
                    selection_reasoning="sandbox 验证通过",
                    confidence="high",
                ),
            )
            if result.success:
                for t in result.new_terminal:
                    print(f"    ✓ {t[:50]}")
                for p in result.new_pending:
                    print(f"    … {p[:50]}")
            else:
                print(f"    ✗ commit failed: {result.error}")
                orch.skip_current(result.error or "commit failed")
        else:
            # 断键失败，尝试 FGI
            fgi_info = orch.explore_fgi()
            if fgi_info.get("n_fgi", 0) > 0:
                fgi_sb = orch.try_fgi(0)
                if fgi_sb.success:
                    print(f"  FGI sandbox: {fgi_sb.template_name}")
                    orch.commit_decision(
                        fgi_template_id=fgi_sb.template_id,
                        fgi_template_name=fgi_sb.template_name,
                    )
                else:
                    orch.skip_current("all sandbox failed")
            else:
                orch.skip_current("sandbox failed, no FGI")

    report = orch.finalize("sandbox-assisted run")
    status = orch.get_status()
    print(f"\n完成: steps={status['steps_executed']}  "
          f"molecules={status['total_molecules']}  "
          f"depth={status['max_depth']}")
    print()
    print(generate_forward_report(orch.get_tree()))
    print("\n✓ Imatinib 沙盒辅助测试通过")


if __name__ == "__main__":
    test_layered_context()
    test_explore_and_sandbox()
    test_explore_fgi()
    test_full_imatinib_with_sandbox()
    print("\n" + "=" * 60)
    print("✓ 所有沙盒测试通过")
