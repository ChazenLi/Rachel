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

from .report_projection import generate_forward_report
from .retro_tree import (
    RetrosynthesisTree,
    ReactionNode,
    MoleculeNode,
    MoleculeRole,
)

# ─────────────────────────────────────────────────────────────────────────
# 起始原料清单
# ─────────────────────────────────────────────────────────────────────────

def _node_id_sort_key(node_id: str) -> tuple[int, str]:
    if not isinstance(node_id, str):
        return (10**9, str(node_id))
    if "_" not in node_id:
        return (10**9, node_id)
    try:
        return (int(node_id.split("_")[-1]), node_id)
    except ValueError:
        return (10**9, node_id)

def get_terminal_list(tree: RetrosynthesisTree) -> List[Dict[str, Any]]:
    """收集所有 terminal 叶节点（起始原料），返回结构化清单。"""
    result: List[Dict[str, Any]] = []
    for node in tree.get_starting_material_nodes():
        entry: Dict[str, Any] = {
            "smiles": node.smiles,
            "node_id": node.node_id,
        }
        if node.complexity:
            entry["cs_score"] = node.cs_score
            entry["classification"] = node.complexity.get("classification", "")
        result.append(entry)
    result.sort(key=lambda x: _node_id_sort_key(str(x.get("node_id", ""))))
    return result


# ─────────────────────────────────────────────────────────────────────────
# 可视化数据
# ─────────────────────────────────────────────────────────────────────────



def _build_visualization_graph(tree: RetrosynthesisTree) -> Dict[str, Any]:
    """Build a tree-shaped display projection from canonical route data.

    The chemistry model deduplicates molecules by canonical SMILES, so a valid
    route can be a DAG when a precursor is reused.  Tree-oriented reports need
    occurrence-specific display nodes while preserving the canonical identity.
    """
    nodes: List[Dict[str, Any]] = []
    edges: List[Dict[str, Any]] = []
    children_map: Dict[str, List[str]] = defaultdict(list)
    reaction_specs: List[Dict[str, Any]] = []
    molecule_sources: Dict[str, str] = {}
    product_to_reactions: Dict[str, List[ReactionNode]] = defaultdict(list)
    for rxn in tree.reaction_nodes:
        product_to_reactions[rxn.product_node].append(rxn)

    used_molecule_ids: set[str] = set()
    used_reaction_ids: set[str] = set()

    def _new_molecule_id(base_id: str, parent_rxn_id: str = "", child_idx: int = 0) -> str:
        if base_id not in used_molecule_ids:
            used_molecule_ids.add(base_id)
            return base_id
        stem = f"{base_id}__{parent_rxn_id}_{child_idx}" if parent_rxn_id else f"{base_id}__occ"
        candidate = stem
        suffix = 2
        while candidate in used_molecule_ids:
            candidate = f"{stem}_{suffix}"
            suffix += 1
        used_molecule_ids.add(candidate)
        return candidate

    def _new_reaction_id(rxn: ReactionNode, product_display_id: str) -> str:
        if rxn.step_id not in used_reaction_ids:
            used_reaction_ids.add(rxn.step_id)
            return rxn.step_id
        stem = f"{rxn.step_id}__{product_display_id}"
        candidate = stem
        suffix = 2
        while candidate in used_reaction_ids:
            candidate = f"{stem}_{suffix}"
            suffix += 1
        used_reaction_ids.add(candidate)
        return candidate

    def _visit(base_id: str, display_id: str, display_depth: int, stack: set[str]) -> None:
        mol = tree.molecule_nodes.get(base_id)
        if mol is None:
            return

        molecule_sources[display_id] = base_id
        node_data: Dict[str, Any] = {
            "id": display_id,
            "canonical_id": base_id,
            "type": "molecule",
            "smiles": mol.smiles,
            "role": mol.role,
            "depth": display_depth,
        }
        if display_id != base_id:
            node_data["display_instance_of"] = base_id
        if mol.complexity:
            node_data["cs_score"] = mol.cs_score
            node_data["classification"] = mol.complexity.get("classification", "")
        nodes.append(node_data)

        # Guard report generation against malformed cyclic route data.
        if base_id in stack:
            return

        for rxn in product_to_reactions.get(base_id, []):
            rxn_display_id = _new_reaction_id(rxn, display_id)
            label = rxn.reaction_type or (
                rxn.template_evidence.template_name
                if rxn.template_evidence else ""
            )
            nodes.append({
                "id": rxn_display_id,
                "canonical_id": rxn.step_id,
                "type": "reaction",
                "label": label,
                "depth": display_depth,
                "reaction_smiles": rxn.reaction_smiles,
            })
            edges.append({
                "source": display_id,
                "target": rxn_display_id,
                "type": "retro_product",
            })

            reactant_display_ids: List[str] = []
            pending_children: List[tuple[str, str]] = []
            for child_idx, rid in enumerate(rxn.reactant_nodes):
                if rid not in tree.molecule_nodes:
                    continue
                child_display_id = _new_molecule_id(rid, rxn_display_id, child_idx)
                reactant_display_ids.append(child_display_id)
                children_map[display_id].append(child_display_id)
                edges.append({
                    "source": rxn_display_id,
                    "target": child_display_id,
                    "type": "retro_reactant",
                })
                pending_children.append((rid, child_display_id))

            for rid, child_display_id in pending_children:
                _visit(rid, child_display_id, display_depth + 1, stack | {base_id})

            reaction_specs.append({
                "id": rxn_display_id,
                "product_id": display_id,
                "reactant_ids": reactant_display_ids,
                "label": label,
            })

    target_node = tree.get_molecule_by_smiles(tree.target)
    root_id = ""
    if target_node is not None:
        root_id = _new_molecule_id(target_node.node_id)
        _visit(target_node.node_id, root_id, target_node.depth, set())

    return {
        "nodes": nodes,
        "edges": edges,
        "root_id": root_id,
        "children_map": dict(children_map),
        "reaction_specs": reaction_specs,
        "molecule_sources": molecule_sources,
    }


def to_visualization_data(tree: RetrosynthesisTree) -> Dict[str, Any]:
    """Output occurrence-expanded nodes/edges for tree-oriented consumers."""
    graph = _build_visualization_graph(tree)
    meta = {
        "target": tree.target,
        "target_name": tree.target_name,
        "status": tree.status,
        "total_steps": tree.total_steps,
        "total_molecules": len(tree.molecule_nodes),
        "display_molecules": sum(1 for n in graph["nodes"] if n["type"] == "molecule"),
    }
    return {"nodes": graph["nodes"], "edges": graph["edges"], "meta": meta}


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
    terminals = list(tree.get_starting_material_nodes())
    terminals.sort(key=lambda n: _node_id_sort_key(n.node_id))
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

