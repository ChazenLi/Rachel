#!/usr/bin/env python
"""快速测试 main 编排系统。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.main import (
    RetrosynthesisOrchestrator,
    generate_forward_report,
    get_terminal_list,
    to_visualization_data,
)

# ── 测试 1: auto_run ──
print("=" * 60)
print("测试 1: auto_run (复杂分子)")
print("=" * 60)

orch = RetrosynthesisOrchestrator(
    "O=C(Nc1cccc2cnccc12)c1cc([N+](=O)[O-])c(Sc2c(Cl)cncc2Cl)s1",
    target_name="Test Amide",
    terminal_cs_threshold=1.3,
    max_depth=4,
)
report = orch.auto_run(verbose=True)

# 正向报告
tree = orch.get_tree()
print()
print(generate_forward_report(tree))

# 起始原料
print()
terminals = get_terminal_list(tree)
print(f"Terminal count: {len(terminals)}")
for t in terminals:
    smi = t["smiles"][:50]
    cs = t.get("cs_score", 0)
    print(f"  {smi}  CS={cs:.2f}")

# 树文本
print()
print(tree.print_tree())

# 可视化数据
vis = to_visualization_data(tree)
print(f"\nVisualization: {len(vis['nodes'])} nodes, {len(vis['edges'])} edges")

# JSON 序列化往返
import json
json_str = tree.to_json()
from Rachel.main.retro_tree import RetrosynthesisTree
tree2 = RetrosynthesisTree.from_json(json_str)
assert tree2.target == tree.target
assert len(tree2.molecule_nodes) == len(tree.molecule_nodes)
assert len(tree2.reaction_nodes) == len(tree.reaction_nodes)
print("\n✓ JSON 序列化往返一致")

# ── 测试 2: 简单分子 (terminal) ──
print("\n" + "=" * 60)
print("测试 2: 简单分子 (应直接 terminal)")
print("=" * 60)

orch2 = RetrosynthesisOrchestrator("c1ccccc1", target_name="Benzene")
report2 = orch2.auto_run(verbose=True)
print(f"Status: {orch2.get_status()}")

# ── 测试 3: LLM 交互模式模拟 ──
print("\n" + "=" * 60)
print("测试 3: LLM 交互模式模拟")
print("=" * 60)

from Rachel.main import LLMDecision

orch3 = RetrosynthesisOrchestrator(
    "CC(=O)Oc1ccccc1C(=O)O",
    target_name="Aspirin",
    terminal_cs_threshold=1.2,
)

ctx = orch3.prepare_next()
print(f"Context: tier={ctx.decision_tier}, smiles={ctx.smiles}")
print(f"  CS={ctx.cs_score:.2f}, is_terminal={ctx.is_terminal}")

if ctx.decision_tier == "standard" and ctx.decision_context:
    bonds = ctx.decision_context.get("disconnectable_bonds", [])
    print(f"  Bonds: {len(bonds)}")
    if bonds:
        b = bonds[0]
        alts = b.get("alternatives", [])
        if alts:
            alt = alts[0]
            print(f"  Best: atoms={b['atoms']}, template={alt['template']}")

            # 模拟 LLM 决策
            decision = LLMDecision(
                selection_reasoning="选择 heuristic_score 最高的断键",
                confidence="high",
            )
            result = orch3.commit_decision(
                bond=tuple(b["atoms"]),
                reaction_type=alt["template"].split("(")[0].strip(),
                template_id=alt.get("template_id", ""),
                template_name=alt["template"],
                llm_decision=decision,
            )
            print(f"  Commit: success={result.success}")
            if result.success:
                print(f"  New terminal: {result.new_terminal}")
                print(f"  New pending: {result.new_pending}")

# 继续处理剩余
while not orch3.is_complete():
    ctx = orch3.prepare_next()
    if ctx is None:
        break
    if ctx.decision_tier == "quick_pass":
        orch3.accept_terminal(reason="auto quick_pass")
    else:
        orch3.skip_current("test: skip remaining")

print(f"\nFinal status: {orch3.get_status()}")
print()
print(generate_forward_report(orch3.get_tree()))

print("\n✓ 所有测试通过")
