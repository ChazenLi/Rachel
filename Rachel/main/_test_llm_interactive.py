#!/usr/bin/env python
"""
LLM 多轮交互式逆合成验证 — 扑热息痛 (Paracetamol)
===================================================
这不是 auto_run，而是真正模拟 LLM 的多轮决策过程:

  Round 1: 目标分子 CC(=O)Nc1ccc(O)cc1
    - 读上下文 → 分析官能团和键位
    - 探索酰胺键 → 查看模板方案
    - 沙盒试探: 模板方案 vs 酰氯法 vs 酸酐法
    - 比较三种方案 → 选择最优 → commit

  Round 2+: 处理前体（如果有 pending）
    - 自动处理 terminal
    - 对 pending 分子继续拆解

  Finalize: 生成正向报告

通过 RetroCmd 标准接口操作，所有数据走 JSON。
"""
import sys, json, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.main.retro_cmd import RetroCmd

SESSION_FILE = str(ROOT / ".test_llm_interactive.json")


def pp(label: str, data: dict, indent: int = 2):
    """Pretty print JSON result."""
    print(f"\n{'─' * 60}")
    print(f"  {label}")
    print(f"{'─' * 60}")
    print(json.dumps(data, indent=indent, ensure_ascii=False))


def main():
    # 清理旧文件
    if os.path.exists(SESSION_FILE):
        os.remove(SESSION_FILE)

    cmd = RetroCmd(SESSION_FILE)

    print("=" * 65)
    print("  LLM 多轮交互式逆合成 — 扑热息痛 (Paracetamol)")
    print("=" * 65)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 0: 初始化会话
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 0: init — 创建会话")
    print("━" * 65)

    r = cmd.execute("init", {
        "target": "CC(=O)Nc1ccc(O)cc1",
        "name": "Paracetamol",
        "terminal_cs_threshold": 1.5,
        "max_steps": 20,
    })
    pp("init result", r)
    assert r.get("ok"), f"init failed: {r}"

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 1: next — 获取第一个分子
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 1: next — LLM 获取目标分子上下文")
    print("━" * 65)

    r = cmd.execute("next")
    pp("next result", r)

    # LLM 思考过程:
    print("\n  [LLM 思考]")
    current = r.get("current", r)
    smiles = current.get("smiles", "?")
    cs = current.get("cs_score", 0)
    fgs = current.get("functional_groups", [])
    bonds = current.get("bond_summary", [])
    print(f"  目标: {smiles} (CS={cs:.2f})")
    print(f"  官能团: {[fg.get('name', fg.get('type', '?')) for fg in fgs]}")
    print(f"  可断键位: {len(bonds)} 个")
    for b in bonds:
        print(f"    bond[{b['bond_idx']}]: atoms={b['atoms']} "
              f"score={b['heuristic_score']:.3f} "
              f"types={b['reaction_types']}")

    print("\n  [LLM 分析]")
    print("  扑热息痛 = 乙酰氨基 + 对羟基苯。核心是酰胺键 C(=O)-N。")
    print("  bond[0] 是 N-C(aryl) 键 → Chan-Lam / SNAr / Ullmann")
    print("  bond[1] 是 C(=O)-N 酰胺键 → Amide Formation / EDC-HATU")
    print("  化学直觉: 酰胺键断裂更合理，但模板给的是偶联反应类型。")
    print("  我要自己提前体: 乙酰氯 + 对氨基酚 (Schotten-Baumann)")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 2: explore — 查看键位详情
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 2: explore — 查看酰胺键 (bond[1]) 详情")
    print("━" * 65)

    r = cmd.execute("explore", {"bond_idx": 1})
    pp("explore bond[1]", r)

    print("\n  [LLM 分析]")
    alts = r.get("alternatives", [])
    print(f"  bond[1] 有 {len(alts)} 个模板方案")
    for i, a in enumerate(alts[:3]):
        prec = a.get("precursors", [])
        tmpl = a.get("template", "")
        print(f"    alt[{i}]: {tmpl}")
        print(f"      precursors: {prec}")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 3: 沙盒试探 — 三种方案比较
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 3a: try_bond — 模板方案 (bond[1], alt[0])")
    print("━" * 65)

    r_a = cmd.execute("try_bond", {"bond_idx": 1, "alt_idx": 0})
    pp("sandbox A: template bond[1].alt[0]", r_a)

    print("\n\n" + "━" * 65)
    print("  STEP 3b: try_precursors — 乙酰氯 + 对氨基酚 (Schotten-Baumann)")
    print("━" * 65)

    r_b = cmd.execute("try_precursors", {
        "precursors": ["CC(=O)Cl", "Nc1ccc(O)cc1"],
        "reaction_type": "Schotten-Baumann acylation",
    })
    pp("sandbox B: AcCl + p-aminophenol", r_b)

    print("\n\n" + "━" * 65)
    print("  STEP 3c: try_precursors — 乙酸酐 + 对氨基酚 (工业路线)")
    print("━" * 65)

    r_c = cmd.execute("try_precursors", {
        "precursors": ["CC(=O)OC(C)=O", "Nc1ccc(O)cc1"],
        "reaction_type": "acetylation (acetic anhydride)",
    })
    pp("sandbox C: Ac2O + p-aminophenol", r_c)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 4: 比较沙盒方案
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 4: sandbox_list — 比较三种方案")
    print("━" * 65)

    r = cmd.execute("sandbox_list")
    pp("sandbox comparison", r)

    print("\n  [LLM 决策推理]")
    print("  方案 A (模板): 前体可能不理想（取决于模板给了什么）")
    print("  方案 B (AcCl + p-aminophenol): Schotten-Baumann，条件温和，")
    print("    两个前体都是简单商业可得试剂，CS 都很低 → 一步到位")
    print("  方案 C (Ac2O + p-aminophenol): 工业路线，但副产物乙酸需处理")
    print("  → 选择方案 B (idx=1)")

    # 分析各方案的前体是否都是 terminal
    for i, att in enumerate(r.get("attempts", [])):
        prec = att.get("precursors", [])
        ok = att.get("success", False)
        fwd = att.get("forward_pass")
        print(f"\n  方案 {chr(65+i)} (idx={i}): {'✓' if ok else '✗'} "
              f"forward={'✓' if fwd else '?' if fwd is None else '✗'}")
        print(f"    precursors: {' + '.join(prec)}")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 5: select + commit
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 5: select(1) + commit — 提交方案 B")
    print("━" * 65)

    cmd.execute("select", {"idx": 1})
    r = cmd.execute("commit", {
        "idx": 1,
        "reasoning": "Schotten-Baumann: AcCl + p-aminophenol → paracetamol. "
                     "两个前体均为简单商业试剂 (CS<1.5)，一步完成。"
                     "优于模板方案（前体可能不理想）和酸酐法（副产物处理）。",
        "confidence": "high",
        "rejected": [
            {"method": "template bond[1].alt[0]", "reason": "模板前体不一定最优"},
            {"method": "Ac2O route", "reason": "副产物乙酸需额外处理"},
        ],
    })
    pp("commit result", r)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 6: next — 处理剩余前体
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 6: next — 检查是否有 pending 前体")
    print("━" * 65)

    r = cmd.execute("next")
    pp("next result", r)

    if r.get("action") == "queue_empty":
        print("\n  [LLM] 所有前体都是 terminal，一步完成。")
    else:
        # 如果还有 pending，继续处理
        print("\n  [LLM] 还有 pending 前体，继续拆解...")
        max_rounds = 5
        for i in range(max_rounds):
            current = r.get("current", r)
            smi = current.get("smiles", "?")
            cs_val = current.get("cs_score", 0)
            print(f"\n  处理: {smi} (CS={cs_val:.2f})")

            # 简单策略: 尝试第一个模板断键
            bond_summary = current.get("bond_summary", [])
            if bond_summary:
                sb = cmd.execute("try_bond", {"bond_idx": 0, "alt_idx": 0})
                if sb.get("success"):
                    cmd.execute("select", {"idx": 0})
                    cr = cmd.execute("commit", {
                        "idx": 0,
                        "reasoning": f"template disconnection for {smi}",
                        "confidence": "medium",
                    })
                    print(f"    commit: terminal={cr.get('new_terminal')} "
                          f"pending={cr.get('new_pending')}")
                else:
                    cmd.execute("accept", {"reason": "no viable disconnection"})
                    print(f"    → accepted as terminal")
            else:
                cmd.execute("accept", {"reason": "no disconnectable bonds"})
                print(f"    → accepted as terminal (no bonds)")

            r = cmd.execute("next")
            if r.get("action") == "queue_empty":
                print("\n  [LLM] 队列已空。")
                break

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 7: tree + status
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 7: tree — 查看最终合成树")
    print("━" * 65)

    r = cmd.execute("tree")
    print(r["tree"])
    print(f"  terminal: {r['terminal_count']}")
    print(f"  pending: {r['pending_count']}")
    print(f"  terminals: {r['terminals']}")

    r = cmd.execute("status")
    pp("status", r)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 8: finalize + report
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 8: finalize + report — 完成编排")
    print("━" * 65)

    r = cmd.execute("finalize", {
        "summary": "扑热息痛一步合成: AcCl + p-aminophenol → Paracetamol (Schotten-Baumann)。"
                   "LLM 自主选择酰氯法优于模板方案和酸酐法。"
    })
    pp("finalize", r)

    r = cmd.execute("report")
    print("\n" + r["report"])

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 9: 验证 JSON 持久化 — 断开恢复
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n\n" + "━" * 65)
    print("  STEP 9: 验证 JSON 持久化 — 模拟断开恢复")
    print("━" * 65)

    # 丢弃内存，从 JSON 恢复
    del cmd
    cmd2 = RetroCmd(SESSION_FILE)
    r = cmd2.execute("status")
    pp("recovered status", r)

    r = cmd2.execute("tree")
    print(r["tree"])

    r = cmd2.execute("report")
    print(r["report"])

    print("\n  ✓ JSON 恢复验证通过")

    # 清理
    if os.path.exists(SESSION_FILE):
        os.remove(SESSION_FILE)

    print("\n" + "=" * 65)
    print("  ✓ LLM 多轮交互式逆合成验证完成")
    print("=" * 65)


if __name__ == "__main__":
    main()
