"""
逆合成结果输出模块
==================
规划完成后，将结果导出到 output/ 目录:

  output/
  └── 20260220_153000_Losartan/
      ├── SYNTHESIS_REPORT.html   # 自包含 HTML 可视化报告（核心输出）
      ├── SYNTHESIS_REPORT.md     # Markdown 报告（带图像引用）
      ├── report.txt              # 正向合成报告（纯文本）
      ├── tree.json               # 完整合成树 JSON
      ├── tree.txt                # 合成树文本渲染
      ├── terminals.json          # 起始原料清单
      ├── visualization.json      # nodes/edges 图数据（供前端）
      ├── session.json            # 完整会话快照（可恢复）
      └── images/                 # 分子/反应/合成树 PNG 图像
          ├── mol_0.png
          ├── rxn_1_reaction.png
          └── synthesis_tree.png
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from .retro_session import RetroSession
from .retro_report import generate_forward_report, get_terminal_list, to_visualization_data


# 项目根目录
ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_ROOT = ROOT / "output"


def _sanitize_name(name: str) -> str:
    """清理名称用于文件夹名。"""
    name = re.sub(r'[<>:"/\\|?*\[\](){}]', '', name)
    name = re.sub(r'\s+', '_', name.strip())
    return name[:40] if name else "unnamed"


def _auto_name_molecule(smiles: str) -> str:
    """根据 SMILES 自动生成简短描述名。"""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "molecule"
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        ha = mol.GetNumHeavyAtoms()
        return f"{formula}_{ha}atoms"
    except Exception:
        return "molecule"


def export_results(
    session: RetroSession,
    output_dir: Optional[str] = None,
    name: Optional[str] = None,
) -> Dict[str, Any]:
    """导出逆合成结果到 output/ 目录。

    生成完整可视化报告：HTML + MD + 图像 + JSON 数据。

    Args:
        session: 已完成（或进行中）的 RetroSession
        output_dir: 自定义输出目录（默认自动生成）
        name: 分子名称（默认从 session 或自动命名）

    Returns:
        {"output_dir": str, "files": [...], "summary": str}
    """
    tree = session.orch.tree
    mol_name = name or tree.target_name or _auto_name_molecule(tree.target)
    safe_name = _sanitize_name(mol_name)

    if output_dir:
        out_path = Path(output_dir)
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = OUTPUT_ROOT / f"{timestamp}_{safe_name}"

    out_path.mkdir(parents=True, exist_ok=True)
    out_str = str(out_path)
    files = []

    # ── 1. 可视化报告（图像 + MD + HTML）──
    vis_result = {}
    try:
        from .retro_visualizer import generate_visual_report
        vis_result = generate_visual_report(tree, out_str, mol_name=mol_name)
        if vis_result.get("md_report"):
            files.append(vis_result["md_report"])
        if vis_result.get("html_report"):
            files.append(vis_result["html_report"])
        if vis_result.get("tree_image"):
            files.append(vis_result["tree_image"])
    except Exception as e:
        vis_result = {"success": False, "error": str(e)}

    # ── 2. 正向合成报告 (txt) ──
    report_text = generate_forward_report(tree)
    report_txt = out_path / "report.txt"
    report_txt.write_text(report_text, encoding="utf-8")
    files.append(str(report_txt))

    # ── 3. 合成树 JSON ──
    tree_json = out_path / "tree.json"
    tree_json.write_text(
        json.dumps(tree.to_dict(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(tree_json))

    # ── 4. 合成树文本 ──
    tree_txt = out_path / "tree.txt"
    tree_txt.write_text(tree.print_tree(), encoding="utf-8")
    files.append(str(tree_txt))

    # ── 5. 起始原料清单 ──
    terminals = get_terminal_list(tree)
    terminals_json = out_path / "terminals.json"
    terminals_json.write_text(
        json.dumps(terminals, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(terminals_json))

    # ── 6. 可视化数据 (nodes/edges) ──
    vis_data = to_visualization_data(tree)
    vis_json = out_path / "visualization.json"
    vis_json.write_text(
        json.dumps(vis_data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(vis_json))

    # ── 7. 完整会话快照 ──
    session_json = out_path / "session.json"
    session_json.write_text(
        json.dumps(session.to_dict(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(session_json))

    # 统计图像数量
    n_mol_imgs = len(vis_result.get("mol_images", {}))
    n_rxn_imgs = len(vis_result.get("rxn_images", {}))
    vis_ok = vis_result.get("success", False)

    summary = (
        f"{mol_name}: {tree.total_steps} 步, "
        f"{len(terminals)} 种起始原料, "
        f"可视化={'✓' if vis_ok else '✗'} "
        f"({n_mol_imgs} 分子图 + {n_rxn_imgs} 反应图), "
        f"输出到 {out_path.relative_to(ROOT)}"
    )

    return {
        "output_dir": str(out_path),
        "relative_dir": str(out_path.relative_to(ROOT)),
        "files": files,
        "n_files": len(files),
        "n_images": n_mol_imgs + n_rxn_imgs + (1 if vis_result.get("tree_image") else 0),
        "html_report": vis_result.get("html_report"),
        "md_report": vis_result.get("md_report"),
        "visualization_ok": vis_ok,
        "visualization_error": vis_result.get("error"),
        "summary": summary,
    }
