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
from typing import Any, Dict, List, Optional, Tuple

from .retro_tree import (
    RetrosynthesisTree,
    MoleculeNode,
    ReactionNode,
    MoleculeRole,
)
from .retro_report import (
    generate_forward_report,
    get_terminal_list,
    to_visualization_data,
    _topological_sort,
)


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


def render_synthesis_tree(
    tree: RetrosynthesisTree,
    output_dir: str,
    mol_images: Optional[Dict[str, str]] = None,
) -> Optional[str]:
    """生成合成树总览图 — 正交直角连线的层级树。

    target 在顶部，每层向下展开前体，直角连线连接。
    每个节点渲染为面板（分子图像 + 标签）。
    """
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return None

    font_label = _load_cjk_font(14)
    font_rxn = _load_cjk_font(12)

    # 构建 parent→children 映射
    children_map: Dict[str, List[str]] = {}
    rxn_for_parent: Dict[str, ReactionNode] = {}
    for rxn in tree.reaction_nodes:
        prod_id = rxn.product_node
        if prod_id not in children_map:
            children_map[prod_id] = []
            rxn_for_parent[prod_id] = rxn
        for rid in rxn.reactant_nodes:
            if rid not in children_map[prod_id]:
                children_map[prod_id].append(rid)

    # 找根节点
    target_node = tree.get_molecule_by_smiles(tree.target)
    if target_node is None:
        return None
    root_id = target_node.node_id

    # 渲染每个分子面板
    panel_w, panel_h = 280, 230
    mol_img_size = (panel_w - 20, panel_h - 50)

    panels: Dict[str, Image.Image] = {}
    for nid, mol_node in tree.molecule_nodes.items():
        panel = Image.new("RGB", (panel_w, panel_h), (255, 255, 255))
        draw = ImageDraw.Draw(panel)

        # 标签行
        cs = mol_node.cs_score
        label = f"{nid} [{mol_node.role.upper()}] CS={cs:.2f}"
        draw.text((8, 4), label, fill=(0, 0, 0), font=font_label)

        # 分子图像
        loaded = False
        if mol_images and nid in mol_images:
            img_path = mol_images[nid]
            if os.path.exists(img_path):
                try:
                    mol_img = Image.open(img_path).convert("RGB")
                    mol_img.thumbnail(mol_img_size)
                    x = (panel_w - mol_img.width) // 2
                    y = 24 + max(0, (panel_h - 28 - mol_img.height) // 2)
                    panel.paste(mol_img, (x, y))
                    loaded = True
                except Exception:
                    pass

        if not loaded:
            try:
                Chem, AllChem, rdMolDraw2D, _ = _ensure_rdkit()
                mol = Chem.MolFromSmiles(mol_node.smiles)
                if mol:
                    AllChem.Compute2DCoords(mol)
                    d = rdMolDraw2D.MolDraw2DCairo(mol_img_size[0], mol_img_size[1])
                    d.DrawMolecule(mol)
                    d.FinishDrawing()
                    png_data = d.GetDrawingText()
                    mol_img = Image.open(io.BytesIO(png_data)).convert("RGB")
                    x = (panel_w - mol_img.width) // 2
                    y = 24 + max(0, (panel_h - 28 - mol_img.height) // 2)
                    panel.paste(mol_img, (x, y))
                    loaded = True
            except Exception:
                pass

        if not loaded:
            smi_short = mol_node.smiles[:30] + "..." if len(mol_node.smiles) > 30 else mol_node.smiles
            draw.text((8, panel_h // 2), smi_short, fill=(80, 80, 80), font=font_label)

        # 边框颜色
        if mol_node.role == MoleculeRole.TARGET.value:
            border_color, border_width = (41, 98, 255), 3
        elif mol_node.role == MoleculeRole.TERMINAL.value:
            border_color, border_width = (56, 142, 60), 2
        else:
            border_color, border_width = (150, 150, 150), 1
        draw.rectangle([0, 0, panel_w - 1, panel_h - 1],
                        outline=border_color, width=border_width)
        panels[nid] = panel

    if not panels:
        return None

    # 递归计算子树宽度
    x_gap, y_gap = 40, 80
    subtree_widths: Dict[str, int] = {}
    _visited: set = set()

    def _subtree_width(nid: str) -> int:
        if nid in subtree_widths:
            return subtree_widths[nid]
        if nid in _visited:
            subtree_widths[nid] = panel_w
            return panel_w
        _visited.add(nid)
        kids = children_map.get(nid, [])
        if not kids:
            subtree_widths[nid] = panel_w
        else:
            total = sum(_subtree_width(k) for k in kids) + x_gap * (len(kids) - 1)
            subtree_widths[nid] = max(panel_w, total)
        return subtree_widths[nid]

    _subtree_width(root_id)

    # 递归布局
    positions: Dict[str, Tuple[int, int]] = {}
    margin = 30
    _visited_layout: set = set()

    def _layout(nid: str, x: int, y: int) -> None:
        if nid in _visited_layout:
            return
        _visited_layout.add(nid)
        sw = subtree_widths.get(nid, panel_w)
        node_x = x + (sw - panel_w) // 2
        positions[nid] = (node_x, y)
        kids = children_map.get(nid, [])
        if not kids:
            return
        child_x = x
        child_y = y + panel_h + y_gap
        for kid in kids:
            _layout(kid, child_x, child_y)
            child_x += subtree_widths.get(kid, panel_w) + x_gap

    _layout(root_id, margin, margin)

    if not positions:
        return None

    max_x = max(px + panel_w for px, _ in positions.values()) + margin
    max_y = max(py + panel_h for _, py in positions.values()) + margin
    canvas = Image.new("RGB", (max(max_x, 400), max(max_y, 300)), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    # 画直角连线
    line_color = (100, 100, 100)
    arrow_color = (60, 60, 60)

    for parent_id, kids in children_map.items():
        if parent_id not in positions:
            continue
        px, py = positions[parent_id]
        pbc = (px + panel_w // 2, py + panel_h)
        child_pos = [(positions[k][0] + panel_w // 2, positions[k][1])
                     for k in kids if k in positions]
        if not child_pos:
            continue
        mid_y = pbc[1] + (child_pos[0][1] - pbc[1]) // 2
        draw.line([pbc, (pbc[0], mid_y)], fill=line_color, width=2)
        if len(child_pos) > 1:
            left_x = min(cx for cx, _ in child_pos)
            right_x = max(cx for cx, _ in child_pos)
            draw.line([(left_x, mid_y), (right_x, mid_y)], fill=line_color, width=2)
        for cx, cy in child_pos:
            if len(child_pos) == 1 and cx != pbc[0]:
                draw.line([(pbc[0], mid_y), (cx, mid_y)], fill=line_color, width=2)
            draw.line([(cx, mid_y), (cx, cy)], fill=line_color, width=2)
            a = 6
            draw.polygon([(cx - a, cy - 1), (cx + a, cy - 1), (cx, cy + a)], fill=arrow_color)

    # 粘贴面板
    for nid, (px, py) in positions.items():
        if nid in panels:
            canvas.paste(panels[nid], (px, py))

    # 标注反应名称
    for rxn in tree.reaction_nodes:
        prod_id = rxn.product_node
        if prod_id not in positions:
            continue
        px, py = positions[prod_id]
        pbc_x = px + panel_w // 2
        pbc_y = py + panel_h
        kid_ids = [r for r in rxn.reactant_nodes if r in positions]
        if not kid_ids:
            continue
        child_y = positions[kid_ids[0]][1]
        mid_y = pbc_y + (child_y - pbc_y) // 2
        rxn_label = rxn.reaction_type or ""
        if rxn.template_evidence and rxn.template_evidence.template_name:
            rxn_label = rxn.template_evidence.template_name
        if rxn_label:
            # 截断过长标签
            if len(rxn_label) > 40:
                rxn_label = rxn_label[:37] + "..."
            draw.text((pbc_x + 5, mid_y - 16), rxn_label,
                      fill=(150, 50, 50), font=font_rxn)

    images_dir = os.path.join(output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)
    tree_path = os.path.join(images_dir, "synthesis_tree.png")
    canvas.save(tree_path)
    return tree_path


# ─────────────────────────────────────────────────────────────────────────
# HTML 报告（自包含，所有图像内嵌 base64）
# ─────────────────────────────────────────────────────────────────────────

def _generate_html_report(
    tree: RetrosynthesisTree,
    mol_name: str,
    tree_image_path: Optional[str] = None,
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

            rxn_type = rxn.reaction_type or ""
            if rxn.template_evidence and rxn.template_evidence.template_name:
                rxn_type = rxn.template_evidence.template_name

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
                                rtext = r if isinstance(r, str) else str(r)
                                h.append(f'<span class="rejected-tag">{_esc(rtext)}</span>')
                        h.append('</div>')
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

    h.append('</div></body></html>')
    return "\n".join(h)


def _esc(s: str) -> str:
    """HTML 转义。"""
    return (s.replace("&", "&amp;").replace("<", "&lt;")
             .replace(">", "&gt;").replace('"', "&quot;"))


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

            rxn_type = rxn.reaction_type or ""
            if rxn.template_evidence and rxn.template_evidence.template_name:
                rxn_type = rxn.template_evidence.template_name

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

            # 反应图像
            if rxn.step_id in rxn_images:
                md.append(f"\n![反应 {rxn.step_id}]({_fix_path(rxn_images[rxn.step_id], report_dir)})\n")

            # 产物
            md.append(f"**产物**: `{product_smi}`\n")
            if rxn.product_node in mol_images:
                md.append(f"![{rxn.product_node}]({_fix_path(mol_images[rxn.product_node], report_dir)})\n")

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
                        rtext = r if isinstance(r, str) else str(r)
                        md.append(f"  - {rtext}")
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
) -> Dict[str, Any]:
    """一键生成完整可视化报告：图像 + Markdown + HTML。

    生成内容：
      1. 每个分子的 PNG 图像 (images/)
      2. 每步反应的 PNG 图像 (images/)
      3. 合成树总览图 (images/synthesis_tree.png)
      4. Markdown 报告 (SYNTHESIS_REPORT.md)
      5. HTML 报告 (SYNTHESIS_REPORT.html) — 自包含，可直接浏览器打开

    Returns:
        {"success": bool, "files": [...], "error": str or None}
    """
    result: Dict[str, Any] = {
        "success": False,
        "output_dir": output_dir,
        "mol_images": {},
        "rxn_images": {},
        "tree_image": None,
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
        tree_img = render_synthesis_tree(tree, output_dir, mol_images)
        result["tree_image"] = tree_img

        # 4. Markdown 报告
        md_path = os.path.join(output_dir, "SYNTHESIS_REPORT.md")
        text_report = generate_forward_report(tree)
        _write_markdown_report(
            tree, mol_images, rxn_images, tree_img, text_report, md_path, name,
        )
        result["md_report"] = md_path

        # 5. HTML 报告
        html_content = _generate_html_report(tree, name, tree_img)
        html_path = os.path.join(output_dir, "SYNTHESIS_REPORT.html")
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        result["html_report"] = html_path

        result["success"] = True

    except Exception as e:
        result["error"] = str(e)

    return result
