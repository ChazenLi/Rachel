#!/usr/bin/env python
"""
测试 JSON 会话持久化 — 模拟 LLM 完整交互流程
=============================================
演示:
  1. 创建会话 → JSON 文件
  2. LLM 读 JSON 获取上下文
  3. 沙盒试探（多轮），结果写入 JSON
  4. 选中方案 → commit → 写入树
  5. 模拟对话断开 → 从 JSON 恢复
  6. 继续编排直到完成
"""
import sys, json, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.main import RetroSession, generate_forward_report

SESSION_FILE = str(ROOT / ".test_session.json")


def cleanup():
    if os.path.exists(SESSION_FILE):
        os.remove(SESSION_FILE)


def test_full_session():
    cleanup()
    print("=" * 65)
    print("测试: JSON 会话 — 扑热息痛 (Paracetamol) 全流程")
    print("=" * 65)

    # ── 1. 创建会话 ──
    print("\n[1] 创建会话")
    session = RetroSession.create(
        "CC(=O)Nc1ccc(O)cc1",
        session_path=SESSION_FILE,
        target_name="Paracetamol",
        terminal_cs_threshold=1.2,
        max_steps=10,
    )
    print(f"  session_id: {session.session_id}")
    print(f"  saved to: {session.path}")

    # 验证 JSON 文件存在且可读
    with open(SESSION_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    print(f"  target: {data['target']['smiles']}")
    print(f"  status: {data['status']['status']}")

    # ── 2. prepare_next → 获取上下文 ──
    print("\n[2] prepare_next — 获取 LLM 上下文")
    ctx = session.prepare_next()
    if ctx and ctx.get("action") == "auto_terminal":
        print(f"  auto terminal: {ctx['smiles']}")
    elif ctx:
        current = ctx.get("current", ctx)
        print(f"  smiles: {current.get('smiles', '?')}")
        print(f"  cs_score: {current.get('cs_score', '?')}")
        bonds = current.get("bond_summary", [])
        print(f"  bond_summary: {len(bonds)} 键位")
        for bs in bonds[:3]:
            print(f"    [{bs['bond_idx']}] atoms={bs['atoms']} "
                  f"score={bs['heuristic_score']:.3f} "
                  f"types={bs['reaction_types']}")

    # ── 3. 沙盒试探 — 多轮 ──
    print("\n[3] 沙盒试探 — 多轮比较")

    # 方案 A: 模板断键
    print("\n  --- 方案 A: 模板断键 (bond[0], alt[0]) ---")
    sb_a = session.sandbox_try(bond_idx=0, alt_idx=0)
    print(f"  success: {sb_a.get('success')}")
    if sb_a.get("precursors"):
        print(f"  precursors: {sb_a['precursors']}")

    # 方案 B: LLM 自提 — 酰氯法
    print("\n  --- 方案 B: LLM 自提 — 乙酰氯 + 对氨基酚 ---")
    sb_b = session.sandbox_try(
        precursors=["CC(=O)Cl", "Nc1ccc(O)cc1"],
        reaction_type="acylation (acyl chloride)",
    )
    print(f"  success: {sb_b.get('success')}")
    if sb_b.get("precursors"):
        print(f"  precursors: {sb_b['precursors']}")
        for pd in sb_b.get("precursor_details", []):
            tag = "term" if pd["is_terminal"] else "pend"
            print(f"    [{tag}] {pd['smiles']}  CS={pd['cs_score']:.2f}")
    fv = sb_b.get("forward_validation", {})
    if fv:
        print(f"  forward: pass={fv.get('pass')} "
              f"feasibility={fv.get('feasibility_score')}")

    # 方案 C: LLM 自提 — 乙酸酐法（工业路线）
    print("\n  --- 方案 C: LLM 自提 — 乙酸酐 + 对氨基酚 ---")
    sb_c = session.sandbox_try(
        precursors=["CC(=O)OC(C)=O", "Nc1ccc(O)cc1"],
        reaction_type="acetylation (anhydride)",
    )
    print(f"  success: {sb_c.get('success')}")
    if sb_c.get("precursors"):
        print(f"  precursors: {sb_c['precursors']}")

    # 验证 JSON 中沙盒状态
    with open(SESSION_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    sandbox = data.get("current", {}).get("sandbox", {})
    print(f"\n  JSON sandbox: {sandbox['n_attempts']} attempts")
    for i, att in enumerate(sandbox["attempts"]):
        src = att.get("source", "?")
        ok = att.get("success", False)
        prec = att.get("precursors", [])
        print(f"    [{i}] {src}: {'✓' if ok else '✗'} → "
              f"{' + '.join(prec) if prec else 'N/A'}")

    # ── 4. 选中方案 B 并 commit ──
    print("\n[4] 选中方案 B (乙酰氯法) → commit")
    session.sandbox_select(1)  # 方案 B
    result = session.commit(
        attempt_idx=1,
        reasoning="乙酰氯+对氨基酚，Schotten-Baumann 条件温和",
        confidence="high",
        rejected=[
            {"method": "template bond break", "reason": "模板方案前体不理想"},
            {"method": "anhydride", "reason": "副产物乙酸需处理"},
        ],
    )
    print(f"  commit success: {result.get('success')}")
    print(f"  new_terminal: {result.get('new_terminal')}")
    print(f"  new_pending: {result.get('new_pending')}")

    # 验证 commit 后 JSON 状态
    with open(SESSION_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
    current_after = data.get("current")
    tree_steps = data.get("tree", {}).get("total_steps", 0)
    if current_after and isinstance(current_after, dict):
        sandbox_after = current_after.get("sandbox", {})
        print(f"  sandbox after commit: {sandbox_after.get('n_attempts', 0)} attempts (should be 0)")
    else:
        print(f"  current after commit: None (已消费，正常)")
    print(f"  tree total_steps: {tree_steps}")

    # ── 5. 模拟对话断开 → 从 JSON 恢复 ──
    print("\n[5] 模拟对话断开 → 从 JSON 恢复")
    del session  # 丢弃内存中的 session
    session2 = RetroSession.load(SESSION_FILE)
    print(f"  session_id: {session2.session_id}")
    status = session2.orch.get_status()
    print(f"  status: {status['status']}")
    print(f"  steps_executed: {status['steps_executed']}")
    print(f"  pending: {status['pending_count']}")
    terminals = session2.orch.tree.get_terminal_smiles()
    print(f"  terminal: {len(terminals)}")

    # 验证树完整性 — 目标分子应该在树中
    target_node = session2.orch.tree.get_molecule_by_smiles("CC(=O)Nc1ccc(O)cc1")
    assert target_node is not None, "恢复后目标分子丢失!"
    print(f"  target node: {target_node.node_id} role={target_node.role}")

    # ── 6. 继续编排直到完成 ──
    print("\n[6] 继续编排 — 处理剩余 pending 分子")
    max_iters = 8
    for i in range(max_iters):
        ctx = session2.prepare_next()
        if ctx is None:
            print(f"  iter {i}: 队列为空，编排完成")
            break

        if ctx.get("action") == "auto_terminal":
            print(f"  iter {i}: auto_terminal → {ctx['smiles']} "
                  f"(CS={ctx['cs_score']:.2f})")
            continue

        # 模拟 LLM: 尝试第一个模板断键
        current = ctx.get("current", ctx)
        smiles = current.get("smiles", "?")
        bonds = current.get("bond_summary", [])
        print(f"  iter {i}: {smiles} — {len(bonds)} bonds")

        if bonds:
            sb = session2.sandbox_try(bond_idx=0, alt_idx=0)
            if sb.get("success"):
                session2.sandbox_select(0)
                r = session2.commit(
                    attempt_idx=0,
                    reasoning=f"auto-select bond[0] for {smiles}",
                    confidence="medium",
                )
                print(f"    commit: {r.get('success')} "
                      f"terminal={r.get('new_terminal')} "
                      f"pending={r.get('new_pending')}")
            else:
                print(f"    bond[0] failed, accepting as terminal")
                session2.accept_terminal(reason="no viable disconnection")
        else:
            print(f"    no bonds, accepting as terminal")
            session2.accept_terminal(reason="no disconnectable bonds")

    # 最终状态
    final_status = session2.orch.get_status()
    final_terminals = session2.orch.tree.get_terminal_smiles()
    print(f"\n  最终状态: {final_status['status']}")
    print(f"  steps: {final_status['steps_executed']}")
    print(f"  terminal: {len(final_terminals)}")
    print(f"  pending: {final_status['pending_count']}")

    # 保存最终状态并验证 JSON roundtrip
    session2.save()
    with open(SESSION_FILE, "r", encoding="utf-8") as f:
        final_data = json.load(f)
    print(f"  JSON tree nodes: {len(final_data['tree'].get('molecules', {}))}")
    print(f"  JSON seen_smiles: {len(final_data.get('seen_smiles', []))}")

    # 再次 load 验证二次恢复
    session3 = RetroSession.load(SESSION_FILE)
    s3_status = session3.orch.get_status()
    s3_terminals = session3.orch.tree.get_terminal_smiles()
    assert s3_status["steps_executed"] == final_status["steps_executed"], \
        f"二次恢复步数不一致: {s3_status['steps_executed']} vs {final_status['steps_executed']}"
    assert len(s3_terminals) == len(final_terminals), \
        f"二次恢复 terminal 数不一致"
    print(f"  二次恢复验证: ✓ (steps={s3_status['steps_executed']}, "
          f"terminal={len(s3_terminals)})")

    # 打印最终合成树
    print("\n" + "=" * 65)
    print("最终合成树:")
    print("=" * 65)
    print(session3.orch.tree.print_tree())

    # 生成正向报告
    report = generate_forward_report(session3.orch.tree)
    print("\n正向合成报告:")
    print(report)

    cleanup()
    print("\n✓ 全部测试通过")


if __name__ == "__main__":
    test_full_session()
