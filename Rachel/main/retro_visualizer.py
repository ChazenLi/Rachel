"""
合成树可视化
============
为 RetrosynthesisTree 生成分子图像、反应图像、合成树总览图和完整可视化报告。

依赖：rdkit（分子/反应渲染）、PIL（图像拼接）

公开函数：
  - render_molecule_images()  — 为树中每个分子生成 PNG
  - render_reaction_images()  — 为每步反应生成 PNG
  - render_synthesis_tree()   — 生成合成树总览图
  - generate_visual_report()  — 一键生成完整可视化报告（MD + HTML + 图像）
"""

from __future__ import annotations

import base64
import io
import os
import textwrap
from typing import Any, Dict, List, Optional, Tuple

from .retro_tree import (
    RetrosynthesisTree,
    MoleculeNode,
    MoleculeRole,
)
from .retro_report import _build_visualization_graph
from .report_projection import forward_reactions, generate_forward_report


# ─────────────────────────────────────────────────────────────────────────
# RDKit 渲染（延迟导入，缺失时优雅降级）
# ─────────────────────────────────────────────────────────────────────────

def _ensure_rdkit():
    """延迟导入 rdkit。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    return Chem, AllChem, rdMolDraw2D, Draw


def _ensure_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def _load_cjk_font(size: int = 14):
    """加载支持中文的字体。"""
    try:
        from PIL import ImageFont
    except ImportError:
        return None
    candidates = [
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/simhei.ttf",
        "C:/Windows/Fonts/simsun.ttc",
        "/System/Library/Fonts/PingFang.ttc",
        "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
    ]
    for p in candidates:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                continue
    return ImageFont.load_default()


def _mol_image(smiles: str, path: str, legend: str = "",
               size: Tuple[int, int] = (400, 400)) -> bool:
    """渲染单个分子为 PNG。"""
    try:
        Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        AllChem.Compute2DCoords(mol)
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.drawOptions().addStereoAnnotation = True
        d.DrawMolecule(mol, legend=legend)
        d.FinishDrawing()
        _ensure_dir(path)
        d.WriteDrawingText(path)
        return True
    except Exception:
        return False


def _mol_png_bytes(smiles: str, size: Tuple[int, int] = (350, 250)) -> Optional[bytes]:
    """渲染分子为 PNG bytes（用于内嵌 HTML）。"""
    try:
        Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        AllChem.Compute2DCoords(mol)
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.drawOptions().addStereoAnnotation = True
        d.DrawMolecule(mol)
        d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:
        return None


def _rxn_image(rxn_smiles: str, path: str,
               size: Tuple[int, int] = (800, 300)) -> bool:
    """渲染反应 SMILES 为 PNG。"""
    try:
        Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            return False
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.DrawReaction(rxn)
        d.FinishDrawing()
        _ensure_dir(path)
        d.WriteDrawingText(path)
        return True
    except Exception:
        return False


def _rxn_png_bytes(rxn_smiles: str, size: Tuple[int, int] = (700, 200)) -> Optional[bytes]:
    """渲染反应为 PNG bytes。"""
    try:
        Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
        rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
        if rxn is None:
            return None
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.DrawReaction(rxn)
        d.FinishDrawing()
        return d.GetDrawingText()
    except Exception:
        return None


def _to_b64(png_bytes: Optional[bytes]) -> str:
    """PNG bytes → base64 data URI。"""
    if not png_bytes:
        return ""
    return "data:image/png;base64," + base64.b64encode(png_bytes).decode()


# ─────────────────────────────────────────────────────────────────────────
# 公开 API: 图像渲染
# ─────────────────────────────────────────────────────────────────────────

def render_molecule_images(
    tree: RetrosynthesisTree,
    output_dir: str,
) -> Dict[str, str]:
    """为树中每个分子生成 PNG 图像。返回 {node_id: image_path}。"""
    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)

    result: Dict[str, str] = {}
    for nid, mol in tree.molecule_nodes.items():
        cs = mol.cs_score
        legend = f"{mol.role.upper()} CS={cs:.1f}"
        path = os.path.join(images_dir, f"{nid}.png")
        if _mol_image(mol.smiles, path, legend=legend):
            result[nid] = path
    return result


def render_reaction_images(
    tree: RetrosynthesisTree,
    output_dir: str,
) -> Dict[str, str]:
    """为每步反应生成 PNG 图像。返回 {step_id: image_path}。"""
    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)

    result: Dict[str, str] = {}
    for rxn in tree.reaction_nodes:
        if not rxn.reaction_smiles or ">>" not in rxn.reaction_smiles:
            continue
        path = os.path.join(images_dir, f"{rxn.step_id}_reaction.png")
        if _rxn_image(rxn.reaction_smiles, path):
            result[rxn.step_id] = path
    return result


_TREE_MOLECULE_WIDTH = 280
_TREE_MOLECULE_HEIGHT = 230
_TREE_REACTION_WIDTH = 260
_TREE_X_GAP = 40
_TREE_Y_GAP = 26
_TREE_MARGIN = 30
_TREE_MIN_WIDTH = 660


def _wrap_tree_label(value: str, width: int = 34) -> List[str]:
    """Wrap a complete route label without truncating its content."""
    text = " ".join(str(value or "Unspecified reaction").split())
    return textwrap.wrap(
        text,
        width=width,
        break_long_words=True,
        break_on_hyphens=False,
    ) or ["Unspecified reaction"]


def _build_synthesis_tree_layout(tree: RetrosynthesisTree) -> Optional[Dict[str, Any]]:
    """Build one vertical geometry model for both PNG and SVG renderers."""
    graph = _build_visualization_graph(tree)
    root_id = graph["root_id"]
    if not root_id:
        return None

    graph_nodes = {node["id"]: node for node in graph["nodes"]}
    forward_steps = {
        reaction.step_id: index
        for index, reaction in enumerate(forward_reactions(tree), 1)
    }
    nodes: Dict[str, Dict[str, Any]] = {}
    children: Dict[str, List[str]] = {}

    for visual_node in graph["nodes"]:
        if visual_node["type"] != "molecule":
            continue
        nodes[visual_node["id"]] = {
            "id": visual_node["id"],
            "canonical_id": visual_node.get("canonical_id", visual_node["id"]),
            "type": "molecule",
            "width": _TREE_MOLECULE_WIDTH,
            "height": _TREE_MOLECULE_HEIGHT,
            "source_id": graph["molecule_sources"][visual_node["id"]],
        }

    for reaction in graph["reaction_specs"]:
        reaction_node = graph_nodes.get(reaction["id"], {})
        canonical_id = reaction_node.get("canonical_id", reaction["id"])
        label = reaction.get("label") or "Unspecified reaction"
        label_lines = _wrap_tree_label(label)
        nodes[reaction["id"]] = {
            "id": reaction["id"],
            "canonical_id": canonical_id,
            "type": "reaction",
            "width": _TREE_REACTION_WIDTH,
            "height": 38 + (17 * len(label_lines)),
            "step_number": forward_steps.get(canonical_id),
            "label": label,
            "label_lines": label_lines,
        }
        children.setdefault(reaction["product_id"], []).append(reaction["id"])
        children[reaction["id"]] = [
            node_id for node_id in reaction["reactant_ids"] if node_id in nodes
        ]

    if root_id not in nodes:
        return None

    subtree_widths: Dict[str, int] = {}
    visiting: set[str] = set()

    def subtree_width(node_id: str) -> int:
        if node_id in subtree_widths:
            return subtree_widths[node_id]
        node = nodes[node_id]
        if node_id in visiting:
            return int(node["width"])
        visiting.add(node_id)
        child_ids = [child for child in children.get(node_id, []) if child in nodes]
        if child_ids:
            child_width = sum(subtree_width(child) for child in child_ids)
            child_width += _TREE_X_GAP * (len(child_ids) - 1)
        else:
            child_width = 0
        visiting.discard(node_id)
        subtree_widths[node_id] = max(int(node["width"]), child_width)
        return subtree_widths[node_id]

    content_width = subtree_width(root_id)
    canvas_width = max(_TREE_MIN_WIDTH, content_width + (2 * _TREE_MARGIN))
    start_x = (canvas_width - content_width) // 2
    positioned: set[str] = set()

    def position(node_id: str, x: int, y: int) -> None:
        if node_id in positioned:
            return
        positioned.add(node_id)
        node = nodes[node_id]
        span = subtree_widths.get(node_id, int(node["width"]))
        node["x"] = x + ((span - int(node["width"])) // 2)
        node["y"] = y
        child_x = x
        child_y = y + int(node["height"]) + _TREE_Y_GAP
        for child_id in children.get(node_id, []):
            if child_id not in nodes:
                continue
            position(child_id, child_x, child_y)
            child_x += subtree_widths.get(child_id, int(nodes[child_id]["width"]))
            child_x += _TREE_X_GAP

    position(root_id, start_x, _TREE_MARGIN)
    visible_nodes = [node for node in nodes.values() if "x" in node]
    canvas_height = max(
        300,
        max(int(node["y"]) + int(node["height"]) for node in visible_nodes)
        + _TREE_MARGIN,
    )

    segments: List[Tuple[Tuple[int, int], Tuple[int, int]]] = []
    arrows: List[Tuple[int, int]] = []
    for parent_id, child_ids in children.items():
        parent = nodes.get(parent_id)
        visible_children = [nodes[child] for child in child_ids if "x" in nodes.get(child, {})]
        if not parent or "x" not in parent or not visible_children:
            continue
        parent_center = int(parent["x"]) + (int(parent["width"]) // 2)
        parent_bottom = int(parent["y"]) + int(parent["height"])
        child_centers = [
            (int(child["x"]) + (int(child["width"]) // 2), int(child["y"]))
            for child in visible_children
        ]
        mid_y = parent_bottom + ((min(y for _, y in child_centers) - parent_bottom) // 2)
        segments.append(((parent_center, parent_bottom), (parent_center, mid_y)))
        left_x = min(x for x, _ in child_centers)
        right_x = max(x for x, _ in child_centers)
        if left_x != right_x:
            segments.append(((left_x, mid_y), (right_x, mid_y)))
        elif left_x != parent_center:
            segments.append(((parent_center, mid_y), (left_x, mid_y)))
        for child_x, child_y in child_centers:
            segments.append(((child_x, mid_y), (child_x, child_y)))
            arrows.append((child_x, child_y))

    return {
        "root_id": root_id,
        "graph": graph,
        "nodes": nodes,
        "children": children,
        "segments": segments,
        "arrows": arrows,
        "width": canvas_width,
        "height": canvas_height,
    }


def _molecule_png_panel(
    tree: RetrosynthesisTree,
    node: Dict[str, Any],
    mol_images: Optional[Dict[str, str]],
) -> Any:
    from PIL import Image, ImageDraw

    panel_w = int(node["width"])
    panel_h = int(node["height"])
    mol_img_size = (panel_w - 20, panel_h - 50)
    panel = Image.new("RGB", (panel_w, panel_h), (255, 255, 255))
    draw = ImageDraw.Draw(panel)
    font_label = _load_cjk_font(14)
    mol_node = tree.molecule_nodes[node["source_id"]]
    label = f"{node['canonical_id']} [{mol_node.role.upper()}] CS={mol_node.cs_score:.2f}"
    draw.text((8, 4), label, fill=(0, 0, 0), font=font_label)

    loaded = False
    image_path = (mol_images or {}).get(node["source_id"])
    if image_path and os.path.exists(image_path):
        try:
            mol_img = Image.open(image_path).convert("RGB")
            mol_img.thumbnail(mol_img_size)
            image_x = (panel_w - mol_img.width) // 2
            image_y = 24 + max(0, (panel_h - 28 - mol_img.height) // 2)
            panel.paste(mol_img, (image_x, image_y))
            loaded = True
        except Exception:
            pass

    if not loaded:
        try:
            Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
            mol = Chem.MolFromSmiles(mol_node.smiles)
            if mol:
                AllChem.Compute2DCoords(mol)
                drawer = rdMolDraw2D.MolDraw2DCairo(*mol_img_size)
                drawer.drawOptions().addStereoAnnotation = True
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                mol_img = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")
                image_x = (panel_w - mol_img.width) // 2
                image_y = 24 + max(0, (panel_h - 28 - mol_img.height) // 2)
                panel.paste(mol_img, (image_x, image_y))
                loaded = True
        except Exception:
            pass

    if not loaded:
        fallback = _wrap_tree_label(mol_node.smiles, width=30)
        for line_index, line in enumerate(fallback):
            draw.text(
                (8, (panel_h // 2) + (line_index * 16)),
                line,
                fill=(80, 80, 80),
                font=font_label,
            )

    if mol_node.role == MoleculeRole.TARGET.value:
        border_color, border_width = (41, 98, 255), 3
    elif mol_node.role == MoleculeRole.TERMINAL.value:
        border_color, border_width = (56, 142, 60), 2
    else:
        border_color, border_width = (150, 150, 150), 1
    draw.rounded_rectangle(
        [0, 0, panel_w - 1, panel_h - 1],
        radius=6,
        outline=border_color,
        width=border_width,
    )
    return panel


def render_synthesis_tree(
    tree: RetrosynthesisTree,
    output_dir: str,
    mol_images: Optional[Dict[str, str]] = None,
    layout: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Render the vertical route overview as a PNG compatibility asset."""
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return None

    layout = layout or _build_synthesis_tree_layout(tree)
    if not layout:
        return None
    canvas = Image.new("RGB", (layout["width"], layout["height"]), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    for start, end in layout["segments"]:
        draw.line([start, end], fill=(100, 100, 100), width=2)
    for x, y in layout["arrows"]:
        draw.polygon([(x - 6, y - 1), (x + 6, y - 1), (x, y + 6)], fill=(60, 60, 60))

    reaction_font = _load_cjk_font(12)
    reaction_step_font = _load_cjk_font(13)
    for node in layout["nodes"].values():
        if "x" not in node:
            continue
        x, y = int(node["x"]), int(node["y"])
        if node["type"] == "molecule":
            canvas.paste(_molecule_png_panel(tree, node, mol_images), (x, y))
            continue
        width, height = int(node["width"]), int(node["height"])
        draw.rounded_rectangle(
            [x, y, x + width, y + height],
            radius=6,
            fill=(255, 248, 246),
            outline=(190, 104, 88),
            width=2,
        )
        step_label = (
            f"Step {node['step_number']}"
            if node.get("step_number") is not None
            else str(node["canonical_id"])
        )
        draw.text((x + 10, y + 7), step_label, fill=(117, 53, 44), font=reaction_step_font)
        for line_index, line in enumerate(node["label_lines"]):
            draw.text(
                (x + 10, y + 27 + (line_index * 17)),
                line,
                fill=(62, 45, 42),
                font=reaction_font,
            )

    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)
    tree_path = os.path.join(images_dir, "synthesis_tree.png")
    canvas.save(tree_path)
    return tree_path


def _molecule_svg(smiles: str, width: int, height: int) -> Optional[str]:
    try:
        Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg[svg.find("<svg"):] if "<svg" in svg else None
    except Exception:
        return None


def render_synthesis_tree_svg(
    tree: RetrosynthesisTree,
    output_dir: str,
    layout: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Render a genuinely vector route overview for the interactive HTML report."""
    layout = layout or _build_synthesis_tree_layout(tree)
    if not layout:
        return None

    width, height = int(layout["width"]), int(layout["height"])
    chunks = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        (
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
            f"viewBox='0 0 {width} {height}' role='img' aria-labelledby='route-title route-desc'>"
        ),
        "<title id='route-title'>Retrosynthesis route overview</title>",
        "<desc id='route-desc'>Target at the top with reaction steps and precursors expanding downward.</desc>",
        "<style>text{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;letter-spacing:0}.molecule-label{font-size:14px;fill:#20231f}.reaction-step{font-size:13px;font-weight:700;fill:#75352c}.reaction-label{font-size:12px;fill:#3e2d2a}</style>",
        f"<rect width='{width}' height='{height}' fill='#ffffff'/>",
        "<g class='route-connectors' fill='none' stroke='#646464' stroke-width='2'>",
    ]
    for (x1, y1), (x2, y2) in layout["segments"]:
        chunks.append(f"<line x1='{x1}' y1='{y1}' x2='{x2}' y2='{y2}'/>")
    chunks.append("</g><g class='route-arrows' fill='#3c3c3c'>")
    for x, y in layout["arrows"]:
        chunks.append(f"<polygon points='{x - 6},{y - 1} {x + 6},{y - 1} {x},{y + 6}'/>")
    chunks.append("</g>")

    for node in layout["nodes"].values():
        if "x" not in node:
            continue
        x, y = int(node["x"]), int(node["y"])
        node_id = _esc(str(node["id"]))
        canonical_id = _esc(str(node["canonical_id"]))
        if node["type"] == "molecule":
            mol_node = tree.molecule_nodes[node["source_id"]]
            role = _esc(str(mol_node.role))
            if mol_node.role == MoleculeRole.TARGET.value:
                border, border_width = "#2962ff", 3
            elif mol_node.role == MoleculeRole.TERMINAL.value:
                border, border_width = "#388e3c", 2
            else:
                border, border_width = "#969696", 1
            chunks.append(
                f"<g data-node-id='{node_id}' data-canonical-id='{canonical_id}' "
                f"data-node-type='molecule' data-role='{role}' data-x='{x}' data-y='{y}'>"
            )
            chunks.append(
                f"<rect x='{x}' y='{y}' width='{node['width']}' height='{node['height']}' "
                f"rx='6' fill='#ffffff' stroke='{border}' stroke-width='{border_width}'/>"
            )
            label = _esc(
                f"{node['canonical_id']} [{str(mol_node.role).upper()}] CS={mol_node.cs_score:.2f}"
            )
            chunks.append(f"<text class='molecule-label' x='{x + 8}' y='{y + 18}'>{label}</text>")
            mol_width = int(node["width"]) - 20
            mol_height = int(node["height"]) - 50
            molecule_svg = _molecule_svg(mol_node.smiles, mol_width, mol_height)
            if molecule_svg:
                molecule_svg = molecule_svg.replace(
                    "<svg ",
                    f"<svg x='{x + 10}' y='{y + 30}' ",
                    1,
                )
                chunks.append(molecule_svg)
            else:
                fallback_lines = _wrap_tree_label(mol_node.smiles, width=30)
                chunks.append(f"<text class='reaction-label' x='{x + 10}' y='{y + 70}'>")
                for line_index, line in enumerate(fallback_lines):
                    dy = 0 if line_index == 0 else 17
                    chunks.append(
                        f"<tspan x='{x + 10}' dy='{dy}'>{_esc(line)}</tspan>"
                    )
                chunks.append("</text>")
            chunks.append("</g>")
            continue

        step_label = (
            f"Step {node['step_number']}"
            if node.get("step_number") is not None
            else str(node["canonical_id"])
        )
        chunks.append(
            f"<g data-node-id='{node_id}' data-canonical-id='{canonical_id}' "
            f"data-node-type='reaction' data-x='{x}' data-y='{y}'>"
        )
        chunks.append(
            f"<rect x='{x}' y='{y}' width='{node['width']}' height='{node['height']}' "
            "rx='6' fill='#fff8f6' stroke='#be6858' stroke-width='2'/>")
        chunks.append(
            f"<text class='reaction-step' x='{x + 10}' y='{y + 18}'>{_esc(step_label)}</text>"
        )
        chunks.append(f"<text class='reaction-label' x='{x + 10}' y='{y + 38}'>")
        for line_index, line in enumerate(node["label_lines"]):
            dy = 0 if line_index == 0 else 17
            chunks.append(f"<tspan x='{x + 10}' dy='{dy}'>{_esc(line)}</tspan>")
        chunks.append("</text></g>")

    chunks.append("</svg>")
    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)
    tree_path = os.path.join(images_dir, "synthesis_tree.svg")
    with open(tree_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(chunks))
    return tree_path


# ─────────────────────────────────────────────────────────────────────────
# HTML 报告（自包含，所有图像内嵌 base64）
# ─────────────────────────────────────────────────────────────────────────

def _scouting_review_rows(
    scouting_view: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    reviews = (scouting_view or {}).get("latest_review_by_node", {}) or {}
    counts = (scouting_view or {}).get("round_count_by_node", {}) or {}
    seeds = (scouting_view or {}).get("deferred_review_seeds", {}) or {}
    adopted_by_node: Dict[str, List[Dict[str, Any]]] = {}
    for item in (scouting_view or {}).get("adopted_sources", []) or []:
        adopted_by_node.setdefault(str(item.get("node_id", "")), []).append(item)
    return [
        {
            "node_id": node_id,
            "round_count": counts.get(node_id, 0),
            "review": review,
            "seeds": seeds.get(node_id, []) or [],
            "adopted": adopted_by_node.get(node_id, []),
        }
        for node_id, review in reviews.items()
    ]


def _append_scouting_html(
    html: List[str],
    scouting_view: Optional[Dict[str, Any]],
) -> None:
    rows = _scouting_review_rows(scouting_view)
    if not rows:
        return
    html.append('<div class="section">')
    html.append('<h2>Scouting Reviews</h2>')
    for row in rows:
        review = row["review"]
        html.append('<div class="decision-audit">')
        html.append(
            '<div class="audit-row"><strong>Node</strong>: '
            f'{_esc(row["node_id"])} | rounds={row["round_count"]} | '
            f'latest={_esc(review.get("round_id", ""))}</div>'
        )
        html.append(
            '<div class="audit-row"><strong>Review</strong>: '
            f'{_esc(review.get("review_summary", ""))}</div>'
        )
        for item in review.get("shortlisted_directions", []) or []:
            direction = item.get("primary_direction", {}) or {}
            html.append(
                '<div class="audit-row"><strong>Shortlist</strong>: '
                f'{_esc(item.get("task_id", ""))} - '
                f'{_esc(direction.get("direction_summary", ""))}</div>'
            )
            disclosure = _frontier_disclosure(direction)
            if disclosure:
                html.append(
                    f'<div class="audit-row"><strong>Frontier hypothesis</strong>: {_esc(disclosure)}</div>'
                )
        for seed in row["seeds"]:
            html.append(
                '<div class="audit-row"><strong>Deferred seed</strong>: '
                f'{_esc(seed.get("seed_id", ""))} - '
                f'{_esc(seed.get("direction_summary", ""))}</div>'
            )
            disclosure = _frontier_disclosure(seed)
            if disclosure:
                html.append(
                    f'<div class="audit-row"><strong>Frontier hypothesis</strong>: {_esc(disclosure)}</div>'
                )
        if row["adopted"]:
            task_ids = [
                (item.get("scouting_source") or {}).get("task_id", "")
                for item in row["adopted"]
            ]
            html.append(
                '<div class="audit-row"><strong>Adopted tasks</strong>: '
                f'{_esc(", ".join(task_ids))}</div>'
            )
            for item in row["adopted"]:
                disclosure = _frontier_disclosure(item)
                if disclosure:
                    task_id = (item.get("scouting_source") or {}).get("task_id", "")
                    html.append(
                        f'<div class="audit-row"><strong>Adopted frontier { _esc(task_id) }</strong>: {_esc(disclosure)}</div>'
                    )
        else:
            html.append(
                '<div class="audit-row">reviewed; no scouting direction adopted</div>'
            )
        html.append('</div>')
    html.append('</div>')


def _append_scouting_markdown(
    markdown: List[str],
    scouting_view: Optional[Dict[str, Any]],
) -> None:
    rows = _scouting_review_rows(scouting_view)
    if not rows:
        return
    markdown.append("## Scouting Reviews\n")
    for row in rows:
        review = row["review"]
        markdown.append(
            f"### Node {row['node_id']} "
            f"({row['round_count']} round(s), latest {review.get('round_id', '')})\n"
        )
        markdown.append(f"- **Review**: {review.get('review_summary', '')}")
        for item in review.get("shortlisted_directions", []) or []:
            direction = item.get("primary_direction", {}) or {}
            markdown.append(
                f"- **Shortlist {item.get('task_id', '')}**: "
                f"{direction.get('direction_summary', '')}"
            )
            disclosure = _frontier_disclosure(direction)
            if disclosure:
                markdown.append(f"- **Frontier hypothesis**: {disclosure}")
        for seed in row["seeds"]:
            markdown.append(
                f"- **Deferred seed {seed.get('seed_id', '')}**: "
                f"{seed.get('direction_summary', '')}"
            )
            disclosure = _frontier_disclosure(seed)
            if disclosure:
                markdown.append(f"- **Frontier hypothesis**: {disclosure}")
        if row["adopted"]:
            task_ids = [
                (item.get("scouting_source") or {}).get("task_id", "")
                for item in row["adopted"]
            ]
            markdown.append(f"- **Adopted tasks**: {', '.join(task_ids)}")
            for item in row["adopted"]:
                disclosure = _frontier_disclosure(item)
                if disclosure:
                    task_id = (item.get("scouting_source") or {}).get("task_id", "")
                    markdown.append(f"- **Adopted frontier {task_id}**: {disclosure}")
        else:
            markdown.append("- reviewed; no scouting direction adopted")
        markdown.append("")


def _generate_html_report(
    tree: RetrosynthesisTree,
    mol_name: str,
    tree_image_path: Optional[str] = None,
    scouting_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
) -> str:
    """生成自包含 HTML 报告，所有分子/反应图像内嵌为 base64。"""

    sorted_rxns = _topological_sort(tree)
    terminals = get_terminal_list(tree)

    # 预渲染所有分子图像为 base64
    mol_b64: Dict[str, str] = {}
    for nid, mol in tree.molecule_nodes.items():
        png = _mol_png_bytes(mol.smiles, size=(300, 200))
        if png:
            mol_b64[nid] = _to_b64(png)

    # 预渲染所有反应图像
    rxn_b64: Dict[str, str] = {}
    for rxn in tree.reaction_nodes:
        if rxn.reaction_smiles and ">>" in rxn.reaction_smiles:
            png = _rxn_png_bytes(rxn.reaction_smiles, size=(650, 180))
            if png:
                rxn_b64[rxn.step_id] = _to_b64(png)

    # 合成树图像
    tree_b64 = ""
    if tree_image_path and os.path.exists(tree_image_path):
        with open(tree_image_path, "rb") as f:
            tree_b64 = "data:image/png;base64," + base64.b64encode(f.read()).decode()

    # 构建 HTML
    h: List[str] = []
    h.append(_HTML_HEAD.replace("{{TITLE}}", f"{mol_name} — 逆合成报告"))

    # 头部信息
    target_b64 = ""
    target_node = tree.get_molecule_by_smiles(tree.target)
    if target_node and target_node.node_id in mol_b64:
        target_b64 = mol_b64[target_node.node_id]

    h.append('<div class="header">')
    h.append(f'<h1>{_esc(mol_name)} — 逆合成规划报告</h1>')
    h.append(f'<div class="meta">目标: <code>{_esc(tree.target)}</code> | '
             f'状态: <span class="badge badge-{tree.status}">{tree.status}</span> | '
             f'步数: {tree.total_steps} | 深度: {tree.max_depth} | '
             f'起始原料: {len(terminals)} 种</div>')
    if target_b64:
        h.append(f'<img src="{target_b64}" class="target-img" alt="目标分子">')
    h.append('</div>')
    _append_scouting_html(h, scouting_view)

    # 合成树总览
    if tree_b64:
        h.append('<div class="section">')
        h.append('<h2>🌳 合成路线总览</h2>')
        h.append(f'<div class="tree-container"><img src="{tree_b64}" alt="合成树"></div>')
        h.append('</div>')

    # 文本树
    h.append('<div class="section">')
    h.append('<h2>📋 合成树（文本）</h2>')
    h.append(f'<pre class="tree-text">{_esc(tree.print_tree())}</pre>')
    h.append('</div>')

    # 起始原料
    if terminals:
        h.append('<div class="section">')
        h.append(f'<h2>🧪 起始原料 ({len(terminals)} 种)</h2>')
        h.append('<div class="card-grid">')
        for t in terminals:
            smi = t["smiles"]
            nid = t["node_id"]
            cs = t.get("cs_score", 0)
            cls_ = t.get("classification", "")
            img = mol_b64.get(nid, "")
            h.append('<div class="mol-card terminal">')
            if img:
                h.append(f'<img src="{img}" alt="{_esc(smi)}">')
            h.append(f'<div class="mol-info">')
            h.append(f'<code>{_esc(smi)}</code>')
            h.append(f'<span class="cs">CS={cs:.2f} [{cls_}]</span>')
            h.append(f'</div></div>')
        h.append('</div></div>')

    # 正向合成步骤
    if sorted_rxns:
        h.append('<div class="section">')
        h.append(f'<h2>⚗️ 正向合成步骤 ({len(sorted_rxns)} 步)</h2>')

        for i, rxn in enumerate(sorted_rxns, 1):
            product_node = tree.molecule_nodes.get(rxn.product_node)
            product_smi = product_node.smiles if product_node else rxn.product_node

            rxn_type = rxn.reaction_type or (
                rxn.template_evidence.template_name
                if rxn.template_evidence else ""
            )

            h.append(f'<div class="step-card">')
            h.append(f'<div class="step-header">Step {i}: {_esc(rxn_type)}</div>')

            # 反应图像
            if rxn.step_id in rxn_b64:
                h.append(f'<div class="rxn-img"><img src="{rxn_b64[rxn.step_id]}" alt="反应"></div>')

            # 前体
            h.append('<div class="reactants">')
            h.append('<span class="label">前体:</span>')
            for j, rid in enumerate(rxn.reactant_nodes):
                rnode = tree.molecule_nodes.get(rid)
                if rnode:
                    rimg = mol_b64.get(rid, "")
                    h.append(f'<div class="mini-mol">')
                    if rimg:
                        h.append(f'<img src="{rimg}" alt="{_esc(rnode.smiles)}">')
                    h.append(f'<code>{_esc(rnode.smiles)}</code>')
                    h.append(f'<span class="role-tag {rnode.role}">{rnode.role}</span>')
                    h.append(f'</div>')
                    if j < len(rxn.reactant_nodes) - 1:
                        h.append('<span class="plus">+</span>')
            h.append('</div>')

            # 产物
            h.append('<div class="product">')
            h.append('<span class="label">产物:</span>')
            prod_img = mol_b64.get(rxn.product_node, "")
            h.append(f'<div class="mini-mol">')
            if prod_img:
                h.append(f'<img src="{prod_img}" alt="{_esc(product_smi)}">')
            h.append(f'<code>{_esc(product_smi)}</code>')
            h.append(f'</div></div>')

            # 安全声明与结构化工艺条件覆盖
            risk_text = (
                rxn.llm_decision.risk_assessment
                if rxn.llm_decision and rxn.llm_decision.risk_assessment
                else "未显式评估；执行前必须核对精确条件、SDS、SOP 与规模化风险"
            )
            process_conditions = (
                getattr(rxn.llm_decision, "process_conditions", {})
                if rxn.llm_decision else {}
            )
            process_coverage = assess_process_conditions(process_conditions)
            h.append('<div class="safety">')
            h.append(f'<div><strong>安全性说明</strong>: {_esc(risk_text)}</div>')
            h.append(
                f'<div><strong>工艺条件声明</strong>: '
                f'{_esc(format_process_conditions(process_conditions))}</div>'
            )
            h.append(
                f'<div><strong>工艺条件覆盖</strong>: '
                f'{_esc(process_coverage["process_assessment_status"])} '
                f'({_esc(format_process_condition_gaps(process_coverage))})</div>'
            )
            h.append('</div>')

            # 决策信息
            if rxn.llm_decision:
                reasoning = rxn.llm_decision.selection_reasoning or ""
                confidence = rxn.llm_decision.confidence or ""
                rejected = rxn.llm_decision.rejected_alternatives or []
                if reasoning:
                    h.append(f'<div class="decision">')
                    h.append(f'<div class="reasoning">💡 {_esc(reasoning)}</div>')
                    if confidence:
                        h.append(f'<span class="badge badge-{confidence}">{confidence}</span>')
                    if rejected:
                        h.append('<div class="rejected">被拒绝: ')
                        if isinstance(rejected, str):
                            h.append(f'<span class="rejected-tag">{_esc(rejected)}</span>')
                        else:
                            for r in rejected:
                                rtext = _format_rejected_alt(r)
                                h.append(f'<span class="rejected-tag">{_esc(rtext)}</span>')
                        h.append('</div>')
                    h.extend(_decision_audit_html(rxn.llm_decision))
                    h.append('</div>')

            h.append('</div>')  # step-card

        h.append('</div>')  # section

    # 分子一览表
    h.append('<div class="section">')
    h.append('<h2>📊 分子一览</h2>')
    h.append('<table class="mol-table"><thead><tr>')
    h.append('<th>ID</th><th>结构</th><th>SMILES</th><th>角色</th><th>CS</th>')
    h.append('</tr></thead><tbody>')
    for nid, mol in tree.molecule_nodes.items():
        img = mol_b64.get(nid, "")
        cs = mol.cs_score
        h.append('<tr>')
        h.append(f'<td>{nid}</td>')
        h.append(f'<td class="img-cell">{"<img src=" + chr(34) + img + chr(34) + ">" if img else "-"}</td>')
        smi_display = mol.smiles if len(mol.smiles) <= 40 else mol.smiles[:37] + "..."
        h.append(f'<td><code>{_esc(smi_display)}</code></td>')
        h.append(f'<td><span class="role-tag {mol.role}">{mol.role}</span></td>')
        h.append(f'<td>{cs:.2f}</td>')
        h.append('</tr>')
    h.append('</tbody></table></div>')

    # LLM 总结
    if tree.llm_summary:
        h.append('<div class="section">')
        h.append('<h2>🤖 LLM 总结</h2>')
        h.append(f'<p>{_esc(tree.llm_summary)}</p>')
        h.append('</div>')

    if post_route_audit:
        from Rachel.tools.post_route_audit import render_post_route_audit_html

        h.append(render_post_route_audit_html(post_route_audit))

    h.append('</div></body></html>')
    return "\n".join(h)


def _esc(s: str) -> str:
    """HTML 转义。"""
    return (s.replace("&", "&amp;").replace("<", "&lt;")
             .replace(">", "&gt;").replace('"', "&quot;"))


def _format_rejected_alt(item) -> str:
    if isinstance(item, str):
        return item
    if isinstance(item, dict):
        method = (
            item.get("method")
            or item.get("action_id")
            or item.get("candidate_id")
            or item.get("strategy_id")
            or item.get("reaction_type")
            or item.get("name")
            or "alternative"
        )
        reason = item.get("reason") or item.get("rationale") or item.get("note") or ""
        return f"{method}: {reason}" if reason else str(method)
    return str(item)


def _decision_audit_html(decision) -> List[str]:
    audit = getattr(decision, "decision_audit", {}) or {}
    if not audit:
        return []

    selected_candidate = audit.get("selected_action_id", "") or audit.get("selected_candidate_id", "")
    lines = ['<div class="decision-audit"><strong>Decision Audit</strong>']
    if selected_candidate:
        lines.append(
            '<div class="audit-row">'
            f'Selected action: <code>{_esc(selected_candidate)}</code>'
            '</div>'
        )

    decision_source = audit.get("decision_source", "")
    if decision_source:
        lines.append(
            '<div class="audit-row">'
            f'<strong>Decision source</strong>: {_esc(decision_source)}'
            '</div>'
        )

    scouting_source = audit.get("scouting_source") or {}
    if scouting_source:
        source_text = (
            f"round={scouting_source.get('round_id', '')}, "
            f"task={scouting_source.get('task_id', '')}"
        )
        if scouting_source.get("adoption_reason"):
            source_text += f", reason={scouting_source['adoption_reason']}"
        lines.append(
            '<div class="audit-row">'
            f'<strong>Scouting source</strong>: {_esc(source_text)}'
            '</div>'
        )

    reagents = list(audit.get("reagents", []) or [])
    if reagents:
        lines.append(
            '<div class="audit-row">'
            f'<strong>Current-step reagents</strong>: {_esc(", ".join(reagents))}'
            '</div>'
        )

    route_plan = audit.get("route_plan_snapshot") or audit.get("route_plan_brief") or {}
    if route_plan:
        plan_id = route_plan.get("id", "") or audit.get("route_plan_id", "")
        revision = route_plan.get("revision", audit.get("route_plan_revision", ""))
        thesis = route_plan.get("route_thesis", "")
        label = f"{plan_id} r{revision}" if plan_id != "" and revision != "" else plan_id
        route_text = f"{label}: {thesis}" if label else thesis
        alignment = audit.get("route_plan_alignment", "")
        if alignment:
            route_text += f" [{alignment}]"
        if audit.get("route_plan_note"):
            route_text += f" | {audit.get('route_plan_note')}"
        lines.append(
            '<div class="audit-row">'
            f'<strong>Route plan</strong>: {_esc(route_text)}'
            '</div>'
        )

    guidance = audit.get("chemist_guidance_summary", []) or []
    if guidance:
        lines.append('<div class="audit-row"><strong>Chemist guidance</strong>: ')
        for item in guidance[:2]:
            label = item.get("id", "")
            intent = item.get("intent", "")
            if intent:
                label = f"{label} ({intent})"
            lines.append(
                f'<span class="audit-candidate">{_esc(label)}: {_esc(item.get("summary", ""))}</span>'
            )
        lines.append('</div>')

    route_strategy = audit.get("route_strategy_brief") or {}
    if route_strategy:
        sketch_id = route_strategy.get("id", "") or audit.get("route_sketch_id", "")
        macro = route_strategy.get("macro_strategy", "")
        next_step = route_strategy.get("next_executable_step", "")
        route_text = f"{sketch_id}: {macro}" if sketch_id else macro
        if next_step:
            route_text += f" -> {next_step}"
        lines.append(
            '<div class="audit-row">'
            f'<strong>Route sketch</strong>: {_esc(route_text)}'
            '</div>'
        )

    gate = audit.get("validation_gate") or {}
    if gate:
        gate_state = _gate_state({"validation_gate": gate})
        lines.append(
            '<div class="audit-row">'
            f'<strong>Validation gate</strong>: {_esc(gate_state + _gate_summary(gate))}'
            '</div>'
        )

    card_ids = audit.get("applied_experience_card_ids", []) or []
    if card_ids:
        lines.append(
            '<div class="audit-row">'
            f'<strong>Applied experience cards</strong>: {_esc(", ".join(card_ids))}'
            '</div>'
        )

    prompt_events = ((audit.get("prompt_state") or {}).get("events", []) or [])
    if prompt_events:
        lines.append(
            '<div class="audit-row">'
            f'<strong>Prompt events</strong>: {_esc(", ".join(prompt_events[:8]))}'
            '</div>'
        )

    selected = audit.get("selected_attempt") or {}
    why_rejected = (
        selected.get("why_existing_actions_rejected", "")
        or selected.get("why_existing_candidates_rejected", "")
    )
    rationale = selected.get("rationale_summary", "")
    if why_rejected or rationale:
        lines.append('<div class="audit-row"><strong>Custom provenance</strong>: ')
        if why_rejected:
            lines.append(f'<span class="audit-candidate">why actions rejected: {_esc(why_rejected)}</span>')
        if rationale:
            lines.append(f'<span class="audit-candidate">rationale: {_esc(rationale)}</span>')
        lines.append('</div>')

    comparison = audit.get("action_comparison", []) or audit.get("candidate_comparison", []) or []
    if comparison:
        lines.append('<div class="audit-row"><strong>Sandbox evidence</strong>: ')
        for item in comparison[:6]:
            lines.append(
                f'<span class="audit-candidate">{_esc(_candidate_audit_line(item))}</span>'
            )
        if len(comparison) > 6:
            lines.append(f'<span class="audit-candidate">+{len(comparison) - 6} more</span>')
        lines.append('</div>')

    rejected = audit.get("rejected_actions", []) or audit.get("rejected_candidates", []) or []
    if rejected:
        lines.append('<div class="audit-row"><strong>Rejected alternatives</strong>: ')
        for item in rejected:
            lines.append(f'<span class="audit-candidate">{_esc(_format_rejected_alt(item))}</span>')
        lines.append('</div>')

    lines.append('</div>')
    return lines


def _decision_audit_markdown(decision) -> List[str]:
    return format_decision_audit_markdown(decision)


_HTML_HEAD = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{{TITLE}}</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background: #f5f7fa; color: #2c3e50; line-height: 1.6;
  }
  .container { max-width: 1200px; margin: 0 auto; padding: 20px; }
  .header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white; padding: 30px; border-radius: 12px; margin-bottom: 24px;
    text-align: center;
  }
  .header h1 { font-size: 1.8em; margin-bottom: 10px; }
  .header .meta { font-size: 0.95em; opacity: 0.9; }
  .header code { background: rgba(255,255,255,0.2); padding: 2px 6px; border-radius: 4px; font-size: 0.85em; }
  .target-img { max-width: 350px; margin-top: 16px; border-radius: 8px; background: white; padding: 8px; }
  .section {
    background: white; border-radius: 10px; padding: 24px;
    margin-bottom: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.06);
  }
  .section h2 { font-size: 1.3em; margin-bottom: 16px; color: #34495e; border-bottom: 2px solid #eee; padding-bottom: 8px; }
  .tree-container { overflow-x: auto; text-align: center; }
  .tree-container img { max-width: 100%; height: auto; }
  .tree-text { background: #f8f9fa; padding: 16px; border-radius: 8px; font-size: 0.85em; overflow-x: auto; line-height: 1.5; }
  .card-grid { display: flex; flex-wrap: wrap; gap: 16px; }
  .mol-card {
    border: 2px solid #e0e0e0; border-radius: 10px; padding: 12px;
    text-align: center; width: 220px; transition: transform 0.2s;
  }
  .mol-card:hover { transform: translateY(-3px); box-shadow: 0 4px 12px rgba(0,0,0,0.1); }
  .mol-card.terminal { border-color: #4caf50; }
  .mol-card img { max-width: 200px; max-height: 140px; }
  .mol-info { margin-top: 8px; }
  .mol-info code { font-size: 0.75em; word-break: break-all; display: block; color: #555; }
  .mol-info .cs { display: block; font-size: 0.8em; color: #888; margin-top: 4px; }
  .step-card {
    border: 1px solid #e0e0e0; border-radius: 10px; padding: 20px;
    margin-bottom: 16px; background: #fafbfc;
  }
  .step-header {
    font-size: 1.1em; font-weight: 600; color: #2962ff;
    margin-bottom: 12px; padding-bottom: 8px; border-bottom: 1px solid #eee;
  }
  .rxn-img { text-align: center; margin: 12px 0; }
  .rxn-img img { max-width: 100%; max-height: 200px; border-radius: 6px; background: white; padding: 4px; }
  .reactants, .product { margin: 8px 0; display: flex; align-items: center; flex-wrap: wrap; gap: 8px; }
  .label { font-weight: 600; color: #555; min-width: 50px; }
  .mini-mol { display: inline-flex; flex-direction: column; align-items: center; padding: 6px; border: 1px solid #eee; border-radius: 6px; background: white; }
  .mini-mol img { max-width: 160px; max-height: 100px; }
  .mini-mol code { font-size: 0.7em; color: #666; margin-top: 4px; max-width: 180px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
  .plus { font-size: 1.4em; font-weight: bold; color: #999; }
  .role-tag {
    display: inline-block; padding: 2px 8px; border-radius: 12px;
    font-size: 0.75em; font-weight: 600; text-transform: uppercase;
  }
  .role-tag.target { background: #e3f2fd; color: #1565c0; }
  .role-tag.terminal { background: #e8f5e9; color: #2e7d32; }
  .role-tag.intermediate { background: #fff3e0; color: #e65100; }
  .decision { margin-top: 12px; padding: 10px; background: #f0f4ff; border-radius: 8px; border-left: 3px solid #667eea; }
  .reasoning { font-size: 0.9em; color: #444; margin-bottom: 6px; }
  .badge {
    display: inline-block; padding: 2px 10px; border-radius: 12px;
    font-size: 0.8em; font-weight: 600;
  }
  .badge-complete, .badge-high { background: #e8f5e9; color: #2e7d32; }
  .badge-in_progress, .badge-medium { background: #fff3e0; color: #e65100; }
  .badge-failed, .badge-low { background: #ffebee; color: #c62828; }
  .rejected { margin-top: 6px; font-size: 0.85em; color: #777; }
  .rejected-tag {
    display: inline-block; background: #ffebee; color: #c62828;
    padding: 1px 6px; border-radius: 4px; margin: 2px; font-size: 0.8em;
  }
  .safety { margin-top: 10px; padding: 8px; background: #fff8e1; border-left: 3px solid #f9a825; font-size: 0.86em; color: #4b3d1f; }
  .safety div + div { margin-top: 4px; }
  .decision-audit { margin-top: 8px; font-size: 0.85em; color: #374151; }
  .audit-row { margin-top: 3px; }
  .audit-candidate { display: inline-block; margin-right: 8px; color: #4b5563; }
  .mol-table { width: 100%; border-collapse: collapse; }
  .mol-table th { background: #f5f7fa; padding: 10px; text-align: left; font-size: 0.9em; border-bottom: 2px solid #ddd; }
  .mol-table td { padding: 8px 10px; border-bottom: 1px solid #eee; vertical-align: middle; }
  .mol-table .img-cell img { max-width: 120px; max-height: 80px; }
  .mol-table code { font-size: 0.8em; }
</style>
</head>
<body><div class="container">
"""


# ─────────────────────────────────────────────────────────────────────────
# Markdown 报告（带图像引用）
# ─────────────────────────────────────────────────────────────────────────

def _fix_path(p: str, base_dir: Optional[str] = None) -> str:
    """将路径转为相对路径并统一正斜杠。"""
    if not p:
        return p
    if base_dir:
        try:
            p = os.path.relpath(p, base_dir)
        except ValueError:
            pass
    return p.replace("\\", "/")


def _write_markdown_report(
    tree: RetrosynthesisTree,
    mol_images: Dict[str, str],
    rxn_images: Dict[str, str],
    tree_image: Optional[str],
    text_report: str,
    output_path: str,
    mol_name: str,
    scouting_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
) -> None:
    """生成 Markdown 格式的可视化合成报告。"""
    report_dir = os.path.dirname(os.path.abspath(output_path))
    sorted_rxns = _topological_sort(tree)
    terminals = get_terminal_list(tree)

    md: List[str] = []

    md.append(f"# {mol_name} — 逆合成规划报告\n")
    md.append(f"**目标分子 SMILES**: `{tree.target}`\n")
    md.append(f"**状态**: {tree.status} | **总步数**: {tree.total_steps} | "
              f"**最大深度**: {tree.max_depth}\n")

    # 目标分子图像
    target_node = tree.get_molecule_by_smiles(tree.target)
    if target_node and target_node.node_id in mol_images:
        md.append(f"![目标分子]({_fix_path(mol_images[target_node.node_id], report_dir)})\n")

    # 合成树总览
    if tree_image:
        md.append("## 合成路线总览\n")
        md.append(f"![合成树]({_fix_path(tree_image, report_dir)})\n")

    _append_scouting_markdown(md, scouting_view)

    # 起始原料
    if terminals:
        md.append(f"## 起始原料 ({len(terminals)} 种)\n")
        md.append("| 编号 | SMILES | CS Score | 分类 | 图像 |")
        md.append("|------|--------|----------|------|------|")
        for i, t in enumerate(terminals, 1):
            smi = t["smiles"]
            cs = t.get("cs_score", 0)
            cls_ = t.get("classification", "-")
            nid = t["node_id"]
            img = f"![{nid}]({_fix_path(mol_images[nid], report_dir)})" if nid in mol_images else "-"
            md.append(f"| {i} | `{smi}` | {cs:.2f} | {cls_} | {img} |")
        md.append("")

    # 正向合成步骤
    if sorted_rxns:
        md.append(f"## 正向合成步骤 ({len(sorted_rxns)} 步)\n")
        for i, rxn in enumerate(sorted_rxns, 1):
            product_node = tree.molecule_nodes.get(rxn.product_node)
            product_smi = product_node.smiles if product_node else rxn.product_node

            rxn_type = rxn.reaction_type or (
                rxn.template_evidence.template_name
                if rxn.template_evidence else ""
            )

            md.append(f"### Step {i}: {rxn_type}\n")
            md.append(f"**反应**: `{rxn.reaction_smiles}`\n")

            # 前体
            md.append("**前体**:\n")
            for rid in rxn.reactant_nodes:
                rnode = tree.molecule_nodes.get(rid)
                if rnode:
                    md.append(f"- `{rnode.smiles}` [{rnode.role}]")
                    if rid in mol_images:
                        md.append(f"  ![{rid}]({_fix_path(mol_images[rid], report_dir)})\n")
            if rxn.reagents:
                md.append(f"\n**当前步试剂**: `{' + '.join(rxn.reagents)}`\n")

            # 反应图像
            if rxn.step_id in rxn_images:
                md.append(f"\n![反应 {rxn.step_id}]({_fix_path(rxn_images[rxn.step_id], report_dir)})\n")

            # 产物
            md.append(f"**产物**: `{product_smi}`\n")
            if rxn.product_node in mol_images:
                md.append(f"![{rxn.product_node}]({_fix_path(mol_images[rxn.product_node], report_dir)})\n")

            risk_text = (
                rxn.llm_decision.risk_assessment
                if rxn.llm_decision and rxn.llm_decision.risk_assessment
                else "未显式评估；执行前必须核对精确条件、SDS、SOP 与规模化风险"
            )
            process_conditions = (
                getattr(rxn.llm_decision, "process_conditions", {})
                if rxn.llm_decision else {}
            )
            process_coverage = assess_process_conditions(process_conditions)
            md.append(f"- **安全性说明**: {risk_text}")
            md.append(
                f"- **工艺条件声明**: {format_process_conditions(process_conditions)}"
            )
            md.append(
                f"- **工艺条件覆盖**: {process_coverage['process_assessment_status']} "
                f"({format_process_condition_gaps(process_coverage)})"
            )

            # 决策
            if rxn.llm_decision and rxn.llm_decision.selection_reasoning:
                md.append(f"- **选择理由**: {rxn.llm_decision.selection_reasoning}")
            if rxn.llm_decision and rxn.llm_decision.confidence:
                md.append(f"- **置信度**: {rxn.llm_decision.confidence}")
            rejected = (rxn.llm_decision.rejected_alternatives or []) if rxn.llm_decision else []
            if rejected:
                md.append("- **被拒绝**:")
                if isinstance(rejected, str):
                    md.append(f"  - {rejected}")
                else:
                    for r in rejected:
                        rtext = _format_rejected_alt(r)
                        md.append(f"  - {rtext}")
            md.extend(_decision_audit_markdown(rxn.llm_decision))
            md.append("")

    # 分子一览
    md.append("## 分子一览\n")
    md.append("| ID | SMILES | 角色 | CS Score | 图像 |")
    md.append("|-----|--------|------|----------|------|")
    for nid, mol in tree.molecule_nodes.items():
        smi_d = mol.smiles if len(mol.smiles) <= 30 else mol.smiles[:27] + "..."
        cs = mol.cs_score
        img = f"[查看]({_fix_path(mol_images[nid], report_dir)})" if nid in mol_images else "-"
        md.append(f"| {nid} | `{smi_d}` | {mol.role} | {cs:.2f} | {img} |")
    md.append("")

    if tree.llm_summary:
        md.append("## LLM 总结\n")
        md.append(tree.llm_summary)
        md.append("")

    if post_route_audit:
        from Rachel.tools.post_route_audit import render_post_route_audit_markdown

        md.append("---")
        md.append(
            render_post_route_audit_markdown(
                post_route_audit,
                heading_level=2,
            ).rstrip()
        )
        md.append("")

    _ensure_dir(output_path)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))


# ─────────────────────────────────────────────────────────────────────────
# 一键生成完整可视化报告
# ─────────────────────────────────────────────────────────────────────────

def generate_visual_report(
    tree: RetrosynthesisTree,
    output_dir: str,
    mol_name: str = "",
    scouting_view: Optional[Dict[str, Any]] = None,
    post_route_audit: Optional[Dict[str, Any]] = None,
    report_view: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """一键生成完整可视化报告：图像 + Markdown + HTML。

    生成内容：
      1. 每个分子的 PNG 图像 (images/)
      2. 每步反应的 PNG 图像 (images/)
      3. 合成树总览图 (images/synthesis_tree.png + synthesis_tree.svg)
      4. Markdown 报告 (SYNTHESIS_REPORT.md)
      5. HTML 报告 (SYNTHESIS_REPORT.html) — 复用 images/ 资源

    Returns:
        {"success": bool, "files": [...], "error": str or None}
    """
    result: Dict[str, Any] = {
        "success": False,
        "output_dir": output_dir,
        "mol_images": {},
        "rxn_images": {},
        "tree_image": None,
        "tree_svg": None,
        "md_report": None,
        "html_report": None,
        "error": None,
    }

    try:
        os.makedirs(output_dir, exist_ok=True)
        name = mol_name or tree.target_name or "molecule"

        # 1. 分子图像
        mol_images = render_molecule_images(tree, output_dir)
        result["mol_images"] = mol_images

        # 2. 反应图像
        rxn_images = render_reaction_images(tree, output_dir)
        result["rxn_images"] = rxn_images

        # 3. 合成树总览
        tree_layout = _build_synthesis_tree_layout(tree)
        tree_img = render_synthesis_tree(
            tree, output_dir, mol_images, layout=tree_layout
        )
        tree_svg = render_synthesis_tree_svg(
            tree, output_dir, layout=tree_layout
        )
        result["tree_image"] = tree_img
        result["tree_svg"] = tree_svg

        # 4. Markdown 报告
        md_path = os.path.join(output_dir, "SYNTHESIS_REPORT.md")
        text_report = generate_forward_report(
            tree, scouting_view=scouting_view, report_view=report_view
        )
        _write_markdown_report(
            tree,
            mol_images,
            rxn_images,
            tree_img,
            text_report,
            md_path,
            name,
            scouting_view=scouting_view,
            post_route_audit=post_route_audit,
            report_view=report_view,
        )
        result["md_report"] = md_path

        # 5. HTML 报告
        html_content = _generate_html_report(
            tree,
            name,
            tree_img,
            tree_svg_path=tree_svg,
            scouting_view=scouting_view,
            post_route_audit=post_route_audit,
            report_view=report_view,
        )
        html_path = os.path.join(output_dir, "SYNTHESIS_REPORT.html")
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        result["html_report"] = html_path

        result["success"] = True

    except Exception as e:
        result["error"] = str(e)

    return result


# Keep image generation here; all human-facing content uses one projection.
from .report_projection import (
    generate_html_report as _generate_html_report,
    write_markdown_report as _write_markdown_report,
)
