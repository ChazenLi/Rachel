"""单条反应可视化 + 逆合成断键分析。

用法:
    # 可视化一条反应 (前体,产物)
    python Rachel/tools/visualize_reaction.py "CC(C)(C)OC(=O)NN.O=C(O)c1ccc(Br)cc1Cl" "CC(C)(C)OC(=O)NNC(=O)c1ccc(Br)cc1Cl"

    # 只可视化一个分子
    python Rachel/tools/visualize_reaction.py "CC(C)(C)OC(=O)NNC(=O)c1ccc(Br)cc1Cl"

输出 HTML 文件到当前目录，自动打开浏览器。
"""

from __future__ import annotations

import argparse
import os
import sys
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import base64
import io


def mol_to_svg(smiles: str, width: int = 400, height: int = 250, highlight_atoms=None, legend: str = "") -> str:
    """将 SMILES 渲染为 SVG 字符串"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f'<div style="color:red">无法解析: {smiles}</div>'
    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.legendFontSize = 20
    if highlight_atoms:
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, legend=legend)
    else:
        drawer.DrawMolecule(mol, legend=legend)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def mol_info(smiles: str) -> dict:
    """提取分子基本信息"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "error": "parse failed"}
    from rdkit.Chem import Descriptors
    return {
        "smiles": smiles,
        "canonical": Chem.MolToSmiles(mol, canonical=True),
        "atoms": mol.GetNumHeavyAtoms(),
        "mw": round(Descriptors.MolWt(mol), 2),
        "rings": mol.GetRingInfo().NumRings(),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
    }


def run_retro_analysis(product_smiles: str) -> dict:
    """对产物跑逆合成分析"""
    from Rachel.chem_tools._rdkit_utils import canonical, parse_mol, tanimoto
    from Rachel.chem_tools.bond_break import execute_disconnection
    from Rachel.chem_tools.template_scan import find_disconnectable_bonds

    result = {"product": product_smiles, "bonds": [], "error": None}

    try:
        scan = find_disconnectable_bonds(product_smiles)
    except Exception as e:
        result["error"] = str(e)
        return result

    if not scan.get("ok"):
        result["error"] = scan.get("error", "scan failed")
        return result

    for bi in scan.get("bonds", []):
        tpls = bi.get("templates", [])
        atoms = bi["atoms"]
        bond_entry = {
            "atoms": atoms,
            "bond_type": bi.get("bond_type", "?"),
            "score": bi.get("heuristic_score", 0),
            "templates": [t["name"] for t in tpls[:3]],
            "break_results": [],
        }

        if tpls:
            try:
                br = execute_disconnection(product_smiles, tuple(atoms), tpls[0]["name"])
                bond_entry["break_results"].append({
                    "success": br.success,
                    "precursors": br.precursors if br.success else [],
                    "error": br.error,
                    "template": tpls[0]["name"],
                })
            except Exception as e:
                bond_entry["break_results"].append({
                    "success": False, "precursors": [], "error": str(e),
                    "template": tpls[0]["name"],
                })

        result["bonds"].append(bond_entry)

    return result


def compare_precursors(true_prec_smiles: str, gen_precursors: list) -> list:
    """比较真实前体与生成前体"""
    from Rachel.chem_tools._rdkit_utils import canonical, parse_mol, tanimoto

    true_frags = true_prec_smiles.split(".")
    true_cans = {canonical(f) for f in true_frags if canonical(f)}

    comparisons = []
    for gen_smi in gen_precursors:
        gc = canonical(gen_smi)
        gen_mol = parse_mol(gen_smi)
        best_sim = 0.0
        best_match = None
        exact = False

        if gc and gc in true_cans:
            exact = True

        for tf in true_frags:
            tm = parse_mol(tf)
            if tm and gen_mol:
                try:
                    sim = tanimoto(gen_mol, tm)
                    if sim > best_sim:
                        best_sim = sim
                        best_match = tf
                except Exception:
                    pass

        comparisons.append({
            "generated": gen_smi,
            "exact_match": exact,
            "best_tanimoto": round(best_sim, 4),
            "best_match_true": best_match,
        })

    return comparisons


def generate_html(precursors_smiles: str | None, product_smiles: str,
                  retro: dict | None, comparisons_by_bond: list | None) -> str:
    """生成完整的 HTML 可视化页面"""

    html_parts = ["""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>反应可视化</title>
<style>
  body { font-family: -apple-system, sans-serif; margin: 20px; background: #f5f5f5; }
  .card { background: white; border-radius: 8px; padding: 20px; margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
  .reaction-row { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; justify-content: center; }
  .mol-box { text-align: center; }
  .arrow { font-size: 36px; color: #666; padding: 0 10px; }
  .plus { font-size: 28px; color: #999; }
  .info-table { border-collapse: collapse; margin: 10px 0; }
  .info-table td, .info-table th { padding: 4px 12px; border: 1px solid #ddd; font-size: 13px; }
  .info-table th { background: #f0f0f0; text-align: left; }
  h2 { color: #333; border-bottom: 2px solid #4CAF50; padding-bottom: 5px; }
  h3 { color: #555; }
  .match { color: #4CAF50; font-weight: bold; }
  .no-match { color: #f44336; }
  .partial { color: #FF9800; }
  .bond-tag { background: #e3f2fd; padding: 2px 8px; border-radius: 4px; font-size: 12px; margin: 2px; display: inline-block; }
  .score { font-weight: bold; }
  svg { max-width: 100%; }
</style></head><body>"""]

    # 反应可视化
    html_parts.append('<div class="card"><h2>反应可视化</h2><div class="reaction-row">')

    if precursors_smiles:
        frags = precursors_smiles.split(".")
        for i, frag in enumerate(frags):
            if i > 0:
                html_parts.append('<span class="plus">+</span>')
            svg = mol_to_svg(frag, 300, 200, legend=f"前体 {i+1}")
            info = mol_info(frag)
            html_parts.append(f'<div class="mol-box">{svg}<div style="font-size:11px;color:#888">{info.get("formula","")} MW={info.get("mw","")}</div></div>')

        html_parts.append('<span class="arrow">→</span>')

    svg = mol_to_svg(product_smiles, 350, 250, legend="产物")
    info = mol_info(product_smiles)
    html_parts.append(f'<div class="mol-box">{svg}<div style="font-size:11px;color:#888">{info.get("formula","")} MW={info.get("mw","")}</div></div>')
    html_parts.append('</div></div>')

    # 产物信息
    html_parts.append('<div class="card"><h2>分子信息</h2>')
    html_parts.append('<table class="info-table"><tr><th>属性</th><th>产物</th></tr>')
    for k, v in info.items():
        html_parts.append(f'<tr><td>{k}</td><td>{v}</td></tr>')
    html_parts.append('</table></div>')

    # 逆合成分析
    if retro:
        html_parts.append(f'<div class="card"><h2>逆合成断键分析</h2>')
        html_parts.append(f'<p>发现 {len(retro["bonds"])} 个可断键位</p>')

        for i, bond in enumerate(retro["bonds"]):
            atoms = bond["atoms"]
            score = bond["score"]
            html_parts.append(f'<div class="card" style="background:#fafafa">')
            html_parts.append(f'<h3>键位 {i+1}: 原子 ({atoms[0]}, {atoms[1]}) '
                            f'<span class="bond-tag">{bond["bond_type"]}</span> '
                            f'<span class="score">score={score}</span></h3>')

            if bond["templates"]:
                html_parts.append(f'<p>匹配模板: {", ".join(bond["templates"])}</p>')

            for br in bond["break_results"]:
                if br["success"]:
                    html_parts.append(f'<p class="match">✓ 断键成功 (模板: {br["template"]})</p>')
                    html_parts.append('<div class="reaction-row">')

                    # 高亮断键位的产物
                    svg_prod = mol_to_svg(product_smiles, 300, 200,
                                         highlight_atoms=atoms, legend="断键位")
                    html_parts.append(f'<div class="mol-box">{svg_prod}</div>')
                    html_parts.append('<span class="arrow">⇒</span>')

                    for j, prec in enumerate(br["precursors"]):
                        if j > 0:
                            html_parts.append('<span class="plus">+</span>')
                        svg_p = mol_to_svg(prec, 280, 200, legend=f"生成前体 {j+1}")
                        html_parts.append(f'<div class="mol-box">{svg_p}<div style="font-size:10px;color:#666">{prec[:60]}</div></div>')

                    html_parts.append('</div>')

                    # 与真实前体比较
                    if comparisons_by_bond and i < len(comparisons_by_bond):
                        comps = comparisons_by_bond[i]
                        if comps:
                            html_parts.append('<table class="info-table"><tr><th>生成前体</th><th>精确匹配</th><th>最佳Tanimoto</th><th>最相似真实前体</th></tr>')
                            for c in comps:
                                cls = "match" if c["exact_match"] else ("partial" if c["best_tanimoto"] > 0.7 else "no-match")
                                html_parts.append(f'<tr><td>{c["generated"][:50]}</td>'
                                                f'<td class="{cls}">{"✓" if c["exact_match"] else "✗"}</td>'
                                                f'<td class="{cls}">{c["best_tanimoto"]}</td>'
                                                f'<td>{(c["best_match_true"] or "")[:50]}</td></tr>')
                            html_parts.append('</table>')
                else:
                    html_parts.append(f'<p class="no-match">✗ 断键失败: {br["error"]}</p>')

            html_parts.append('</div>')

        html_parts.append('</div>')

    html_parts.append('</body></html>')
    return "\n".join(html_parts)


def main():
    parser = argparse.ArgumentParser(description="反应可视化")
    parser.add_argument("precursors_or_product", help="前体SMILES (用.分隔) 或产物SMILES")
    parser.add_argument("product", nargs="?", default=None, help="产物SMILES (可选)")
    parser.add_argument("--no-retro", action="store_true", help="不跑逆合成分析")
    parser.add_argument("--output", "-o", default=None, help="输出文件名")
    args = parser.parse_args()

    if args.product:
        precursors_smiles = args.precursors_or_product
        product_smiles = args.product
    else:
        precursors_smiles = None
        product_smiles = args.precursors_or_product

    print(f"产物: {product_smiles}")
    if precursors_smiles:
        print(f"前体: {precursors_smiles}")

    # 逆合成分析
    retro = None
    comparisons_by_bond = None
    if not args.no_retro:
        print("运行逆合成分析...")
        retro = run_retro_analysis(product_smiles)
        print(f"  发现 {len(retro['bonds'])} 个可断键位")

        if precursors_smiles and retro["bonds"]:
            comparisons_by_bond = []
            for bond in retro["bonds"]:
                bond_comps = []
                for br in bond["break_results"]:
                    if br["success"] and br["precursors"]:
                        bond_comps = compare_precursors(precursors_smiles, br["precursors"])
                comparisons_by_bond.append(bond_comps)

    # 生成 HTML
    html = generate_html(precursors_smiles, product_smiles, retro, comparisons_by_bond)

    out_file = args.output or "reaction_viz.html"
    out_path = Path(out_file)
    out_path.write_text(html, encoding="utf-8")
    print(f"输出: {out_path.resolve()}")

    # 打开浏览器
    webbrowser.open(str(out_path.resolve()))


if __name__ == "__main__":
    main()
