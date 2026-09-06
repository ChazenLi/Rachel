"""
逆合成结果输出模块
==================
规划完成后，将结果导出到 output/ 目录:

  output/
  └── 20260220_153000_Losartan/
      ├── SYNTHESIS_REPORT.html   # HTML 可视化报告（引用 images/ 资源）
      ├── SYNTHESIS_REPORT.md     # Markdown 报告（带图像引用）
      ├── report.txt              # 正向合成报告（纯文本）
      ├── tree.json               # 完整合成树 JSON
      ├── tree.txt                # 合成树文本渲染
      ├── terminals.json          # 起始原料清单
      ├── visualization.json      # nodes/edges 图数据（供前端）
      ├── session.json            # 完整会话快照（可恢复）
      └── images/                 # 分子/反应图与路线总览资源
          ├── mol_0.png
          ├── rxn_1_reaction.png
          ├── synthesis_tree.png  # Markdown/兼容回退
          └── synthesis_tree.svg  # HTML 交互总览
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from .retro_session import RetroSession
from .node_scouting import build_report_view
from .retro_report import generate_forward_report, to_visualization_data
from .public_protocol import project_public_payload
from .report_projection import build_route_report_view
from .report_visualization import annotate_visualization


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


def public_action_tree_dict(
    tree: Any,
    *,
    post_route_audit: Optional[Dict[str, Any]] = None,
    post_route_audit_status: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Return the public tree plus export-only post-route audit projections."""
    payload = project_public_payload(tree.to_dict())
    if post_route_audit_status:
        payload["post_route_audit_status"] = project_public_payload(
            post_route_audit_status
        )
    if post_route_audit:
        payload["post_route_audit"] = project_public_payload(post_route_audit)
    return payload


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
    scouting_view = build_report_view(
        session.orch.audit_state.node_scouting_rounds,
        tree,
    )
    mol_name = name or tree.target_name or _auto_name_molecule(tree.target)
    safe_name = _sanitize_name(mol_name)

    if output_dir:
        out_path = Path(output_dir)
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = OUTPUT_ROOT / f"{timestamp}_{safe_name}"

    if (out_path / "session.json").resolve() == session.path.resolve():
        return {"error": "output_dir would overwrite the active session file"}

    out_path.mkdir(parents=True, exist_ok=True)
    out_str = str(out_path)
    files = []
    audit_export_error = ""
    for sidecar_name in ("post_route_audit.json", "post_route_audit.md"):
        sidecar = out_path / sidecar_name
        try:
            if sidecar.is_file():
                sidecar.unlink()
        except OSError as exc:
            audit_export_error = str(exc)
            break
    current_audit = session.current_post_route_audit()
    audit_status = session.post_route_audit_status()
    experimental_outcomes = list(session.orch.audit_state.experimental_outcomes)
    report_view = build_route_report_view(
        tree,
        route_plan=session._route_plan_brief(),
        scouting_view=scouting_view,
        post_route_audit=current_audit,
        experimental_outcomes=experimental_outcomes,
        variant_info=session._variant,
        decision_history=session.orch.audit_state.decision_history,
    )

    # ── 1. 可视化报告（图像 + MD + HTML）──
    vis_result = {}
    try:
        from .retro_visualizer import generate_visual_report
        vis_result = generate_visual_report(
            tree,
            out_str,
            mol_name=mol_name,
            scouting_view=scouting_view,
            post_route_audit=current_audit,
            report_view=report_view,
        )
        if vis_result.get("md_report"):
            files.append(vis_result["md_report"])
        if vis_result.get("html_report"):
            files.append(vis_result["html_report"])
        if vis_result.get("tree_image"):
            files.append(vis_result["tree_image"])
        if vis_result.get("tree_svg"):
            files.append(vis_result["tree_svg"])
    except Exception as e:
        vis_result = {"success": False, "error": str(e)}

    # ── 2. 正向合成报告 (txt) ──
    report_text = generate_forward_report(tree, report_view=report_view)
    report_txt = out_path / "report.txt"
    report_txt.write_text(report_text, encoding="utf-8")
    files.append(str(report_txt))

    # ── 3. 合成树 JSON ──
    tree_json = out_path / "tree.json"
    tree_json.write_text(
        json.dumps(
            public_action_tree_dict(
                tree,
                post_route_audit=current_audit,
                post_route_audit_status=audit_status,
            ),
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    files.append(str(tree_json))

    # ── 4. 合成树文本 ──
    tree_txt = out_path / "tree.txt"
    tree_txt.write_text(tree.print_tree(), encoding="utf-8")
    files.append(str(tree_txt))

    # ── 5. 起始原料清单 ──
    terminals = report_view["starting_materials"]
    terminals_json = out_path / "terminals.json"
    terminals_json.write_text(
        json.dumps(terminals, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(terminals_json))

    # ── 6. 可视化数据 (nodes/edges) ──
    vis_data = annotate_visualization(to_visualization_data(tree), report_view)
    vis_json = out_path / "visualization.json"
    vis_json.write_text(
        json.dumps(vis_data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(vis_json))

    # ── 7. 去路径知识配置 ──
    knowledge_snapshot = session.knowledge_profile.snapshot()
    knowledge_json = out_path / "knowledge_profile.json"
    knowledge_json.write_text(
        json.dumps(knowledge_snapshot, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(knowledge_json))

    # ── 8. 可共享会话快照 ──
    exported_session = session.to_dict()
    exported_session.pop("knowledge_roots", None)
    session_json = out_path / "session.json"
    session_json.write_text(
        json.dumps(exported_session, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files.append(str(session_json))

    # ── 9. Finalized-route audit sidecars ──
    if current_audit and not audit_export_error:
        try:
            from Rachel.tools.post_route_audit import write_post_route_audit_files

            files.extend(write_post_route_audit_files(current_audit, out_path))
        except Exception as exc:
            audit_export_error = str(exc)

    # 统计图像数量
    n_mol_imgs = len(vis_result.get("mol_images", {}))
    n_rxn_imgs = len(vis_result.get("rxn_images", {}))
    vis_ok = vis_result.get("success", False)
    try:
        display_dir = out_path.relative_to(ROOT)
    except ValueError:
        display_dir = out_path

    summary = (
        f"{mol_name}: {tree.total_steps} steps, "
        f"{len(terminals)} starting materials, "
        f"visualization={'ok' if vis_ok else 'failed'} "
        f"({n_mol_imgs} molecule images + {n_rxn_imgs} reaction images), "
        f"exported to {display_dir}"
    )

    return {
        "output_dir": str(out_path),
        "relative_dir": str(display_dir),
        "files": files,
        "n_files": len(files),
        "n_images": (
            n_mol_imgs
            + n_rxn_imgs
            + (1 if vis_result.get("tree_image") else 0)
            + (1 if vis_result.get("tree_svg") else 0)
        ),
        "html_report": vis_result.get("html_report"),
        "md_report": vis_result.get("md_report"),
        "visualization_ok": vis_ok,
        "visualization_error": vis_result.get("error"),
        "knowledge_profile_hash": knowledge_snapshot["digest"],
        "post_route_audit_status": audit_status,
        "post_route_audit_export_error": audit_export_error,
        "summary": summary,
    }
