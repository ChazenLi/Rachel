"""
合成报告生成器
==============
从 RetrosynthesisTree 生成三种输出:
  1. generate_forward_report() — 正向合成报告（拓扑排序反转）
  2. get_terminal_list() — 起始原料清单
  3. to_visualization_data() — nodes/edges 图数据（供前端可视化）
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List

from .retro_tree import (
    RetrosynthesisTree,
    ReactionNode,
    MoleculeNode,
    MoleculeRole,
)


# ─────────────────────────────────────────────────────────────────────────
# 正向合成报告
# ─────────────────────────────────────────────────────────────────────────

def generate_forward_report(tree: RetrosynthesisTree) -> str:
    """将逆合成树反转为正向合成报告文本。

    拓扑排序: 叶节点反应在前，target 反应在后。
    """
    if not tree.reaction_nodes:
        target_label = tree.target_name or tree.target
        return f"正向合成报告: {target_label}\n\n（无反应步骤）"

    sorted_rxns = _topological_sort(tree)
    terminals = _collect_terminals(tree)

    lines: List[str] = []

    # 标题
    target_label = tree.target_name or tree.target
    lines.append(f"正向合成报告: {target_label}")
    lines.append("=" * 50)
    lines.append("")

    # 起始原料
    lines.append(f"起始原料 ({len(terminals)} 种):")
    for t in terminals:
        cs_info = ""
        if t.complexity:
            cs_info = f"  [CS={t.cs_score:.1f}, {t.complexity.get('classification', '')}]"
        lines.append(f"  • {t.smiles}{cs_info}")
    lines.append("")

    # 合成步骤
    lines.append(f"合成步骤 ({len(sorted_rxns)} 步):")
    lines.append("-" * 40)

    for i, rxn in enumerate(sorted_rxns, 1):
        lines.append("")
        lines.append(f"Step {i}: {rxn.reaction_smiles}")

        reactant_smiles = _get_reactant_smiles(tree, rxn)
        product_smiles = _get_product_smiles(tree, rxn)
        lines.append(f"  前体: {' + '.join(reactant_smiles)}")
        lines.append(f"  产物: {product_smiles}")

        if rxn.reaction_type:
            lines.append(f"  反应类型: {rxn.reaction_type}")

        if rxn.template_evidence and rxn.template_evidence.template_name:
            lines.append(f"  模板: {rxn.template_evidence.template_name}")

        if rxn.llm_decision and rxn.llm_decision.selection_reasoning:
            lines.append(f"  理由: {rxn.llm_decision.selection_reasoning}")

        if rxn.forward_validation:
            fv = rxn.forward_validation
            score = fv.get("assessment", {}).get("feasibility_score")
            passed = fv.get("assessment", {}).get("pass")
            if score is not None:
                icon = "✓" if passed else "✗"
                lines.append(f"  验证: {icon} feasibility={score:.3f}")

    # 总结
    lines.append("")
    lines.append("-" * 40)
    lines.append(f"总计: {len(sorted_rxns)} 步, {len(terminals)} 种起始原料")
    if tree.llm_summary:
        lines.append(f"\nLLM 总结: {tree.llm_summary}")

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────
# 起始原料清单
# ─────────────────────────────────────────────────────────────────────────

def get_terminal_list(tree: RetrosynthesisTree) -> List[Dict[str, Any]]:
    """收集所有 terminal 叶节点，返回结构化清单。"""
    result: List[Dict[str, Any]] = []
    for node in tree.molecule_nodes.values():
        if node.role != MoleculeRole.TERMINAL.value:
            continue
        entry: Dict[str, Any] = {
            "smiles": node.smiles,
            "node_id": node.node_id,
        }
        if node.complexity:
            entry["cs_score"] = node.cs_score
            entry["classification"] = node.complexity.get("classification", "")
        result.append(entry)
    result.sort(key=lambda x: x.get("node_id", ""))
    return result


# ─────────────────────────────────────────────────────────────────────────
# 可视化数据
# ─────────────────────────────────────────────────────────────────────────

def to_visualization_data(tree: RetrosynthesisTree) -> Dict[str, Any]:
    """输出 nodes/edges 图数据，供前端可视化。

    边方向遵循逆合成方向（产物 → 反应 → 前体），前端可自行反转。
    """
    nodes: List[Dict[str, Any]] = []
    edges: List[Dict[str, Any]] = []

    # 分子节点
    for nid, mol in tree.molecule_nodes.items():
        node_data: Dict[str, Any] = {
            "id": nid,
            "type": "molecule",
            "smiles": mol.smiles,
            "role": mol.role,
            "depth": mol.depth,
        }
        if mol.complexity:
            node_data["cs_score"] = mol.cs_score
            node_data["classification"] = mol.complexity.get("classification", "")
        nodes.append(node_data)

    # 反应节点 + 边
    for rxn in tree.reaction_nodes:
        label = rxn.reaction_type
        if rxn.template_evidence and rxn.template_evidence.template_name:
            label = rxn.template_evidence.template_name

        rxn_data: Dict[str, Any] = {
            "id": rxn.step_id,
            "type": "reaction",
            "label": label,
            "depth": rxn.depth,
            "reaction_smiles": rxn.reaction_smiles,
        }
        nodes.append(rxn_data)

        # 产物 → 反应节点
        edges.append({
            "source": rxn.product_node,
            "target": rxn.step_id,
            "type": "retro_product",
        })

        # 反应节点 → 前体
        for rid in rxn.reactant_nodes:
            edges.append({
                "source": rxn.step_id,
                "target": rid,
                "type": "retro_reactant",
            })

    meta = {
        "target": tree.target,
        "target_name": tree.target_name,
        "status": tree.status,
        "total_steps": tree.total_steps,
        "total_molecules": len(tree.molecule_nodes),
    }

    return {"nodes": nodes, "edges": edges, "meta": meta}


# ─────────────────────────────────────────────────────────────────────────
# 内部辅助
# ─────────────────────────────────────────────────────────────────────────

def _topological_sort(tree: RetrosynthesisTree) -> List[ReactionNode]:
    """拓扑排序: 叶节点反应在前，target 反应在后（正向合成顺序）。"""
    product_to_rxn: Dict[str, ReactionNode] = {}
    for rxn in tree.reaction_nodes:
        product_to_rxn[rxn.product_node] = rxn

    # 依赖图: rxn_A 依赖 rxn_B = rxn_A 的某个前体是 rxn_B 的产物
    in_degree: Dict[str, int] = {rxn.step_id: 0 for rxn in tree.reaction_nodes}
    forward: Dict[str, List[str]] = defaultdict(list)

    for rxn in tree.reaction_nodes:
        for rid in rxn.reactant_nodes:
            if rid in product_to_rxn:
                dep_rxn = product_to_rxn[rid]
                in_degree[rxn.step_id] += 1
                forward[dep_rxn.step_id].append(rxn.step_id)

    # Kahn's algorithm
    queue = sorted([sid for sid, deg in in_degree.items() if deg == 0])
    result: List[ReactionNode] = []
    rxn_by_id = {rxn.step_id: rxn for rxn in tree.reaction_nodes}

    while queue:
        current = queue.pop(0)
        result.append(rxn_by_id[current])
        for neighbor in sorted(forward[current]):
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    return result


def _collect_terminals(tree: RetrosynthesisTree) -> List[MoleculeNode]:
    terminals = [
        n for n in tree.molecule_nodes.values()
        if n.role == MoleculeRole.TERMINAL.value
    ]
    terminals.sort(key=lambda n: n.node_id)
    return terminals


def _get_reactant_smiles(tree: RetrosynthesisTree, rxn: ReactionNode) -> List[str]:
    result = []
    for rid in rxn.reactant_nodes:
        node = tree.molecule_nodes.get(rid)
        result.append(node.smiles if node else rid)
    return result


def _get_product_smiles(tree: RetrosynthesisTree, rxn: ReactionNode) -> str:
    node = tree.molecule_nodes.get(rxn.product_node)
    return node.smiles if node else rxn.product_node
