"""低 Tanimoto 案例分析 — 找出逆合成预测最差的案例并分析原因。

用法:
    python Rachel/tools/_near_miss.py                         # BTF1, 200条, 显示最差20条
    python Rachel/tools/_near_miss.py --file BTF1 --size 500 --top 30
    python Rachel/tools/_near_miss.py --threshold 0.5         # 只显示 tanimoto < 0.5 的
    python Rachel/tools/_near_miss.py --html                  # 同时输出 HTML 可视化

分析维度:
  - 产物结构特征 (环数、杂原子、分子量)
  - 模板匹配情况 (匹配了哪些模板、是否走了 fallback)
  - 生成前体 vs 真实前体对比
  - 失败原因分类
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from Rachel.chem_tools._rdkit_utils import canonical, parse_mol, tanimoto
from Rachel.chem_tools.bond_break import execute_disconnection, execute_fgi
from Rachel.chem_tools.template_scan import find_disconnectable_bonds, scan_applicable_reactions
from Rachel.tests.data_driven.loader import load_csv

DATA_DIR = ROOT / "Rachel" / "data"
FILE_MAP = {
    "ATF":  DATA_DIR / "USPTO50K" / "datasetATF.csv",
    **{f"BTF{i}": DATA_DIR / "USPTO50K" / f"datasetBTF{i}.csv" for i in range(1, 11)},
    "NEW":  DATA_DIR / "NEW" / "usptonew.csv",
}


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass
class CaseResult:
    """单条记录的逆合成评估结果"""
    row_index: int
    product: str
    true_precursors: str
    best_tanimoto: float
    exact_set_match: bool
    # 产物特征
    n_rings: int = 0
    n_heavy: int = 0
    mw: float = 0.0
    heteroatoms: str = ""
    formula: str = ""
    # 断键信息
    n_bonds_found: int = 0
    n_bonds_succeeded: int = 0
    n_templates_matched: int = 0
    used_fallback: bool = False
    # 最佳断键结果
    best_gen_precursors: List[str] = field(default_factory=list)
    best_template: str = ""
    best_bond: Tuple[int, int] = (0, 0)
    # 所有断键尝试
    all_attempts: List[Dict[str, Any]] = field(default_factory=list)
    # FGI 尝试
    fgi_attempted: int = 0
    fgi_succeeded: int = 0
    best_fgi_tanimoto: float = 0.0
    best_fgi_precursors: List[str] = field(default_factory=list)
    best_fgi_template: str = ""
    # 综合最佳 (断键 or FGI)
    combined_best_tanimoto: float = 0.0
    # 失败分类
    failure_category: str = ""


def mol_features(smiles: str) -> Dict[str, Any]:
    """提取分子结构特征"""
    mol = parse_mol(smiles)
    if mol is None:
        return {}
    return {
        "n_rings": mol.GetRingInfo().NumRings(),
        "n_heavy": mol.GetNumHeavyAtoms(),
        "mw": round(Descriptors.MolWt(mol), 1),
        "heteroatoms": ",".join(sorted(set(
            a.GetSymbol() for a in mol.GetAtoms() if a.GetSymbol() not in ("C", "H")
        ))),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
    }


# ---------------------------------------------------------------------------
# 核心评估
# ---------------------------------------------------------------------------

def evaluate_single(rec) -> CaseResult:
    """对单条记录做完整逆合成评估，返回详细结果"""
    product = rec.product_smiles
    true_prec = rec.precursors_smiles

    case = CaseResult(
        row_index=rec.row_index,
        product=product,
        true_precursors=true_prec,
        best_tanimoto=0.0,
        exact_set_match=False,
    )

    # 产物特征
    feat = mol_features(product)
    if feat:
        case.n_rings = feat["n_rings"]
        case.n_heavy = feat["n_heavy"]
        case.mw = feat["mw"]
        case.heteroatoms = feat["heteroatoms"]
        case.formula = feat["formula"]

    # 真实前体
    true_frags = true_prec.split(".")
    true_set = {canonical(f) for f in true_frags if canonical(f)}
    true_mols = [m for f in true_frags if (m := parse_mol(f))]

    if not true_set or not true_mols:
        case.failure_category = "true_precursor_parse_failed"
        return case

    # 扫描可断键位
    try:
        scan = find_disconnectable_bonds(product)
    except Exception as e:
        case.failure_category = f"scan_exception: {e}"
        return case

    if not scan.get("ok") or not scan.get("bonds"):
        case.failure_category = "no_disconnectable_bonds"
        return case

    bonds = scan["bonds"]
    case.n_bonds_found = len(bonds)

    # 对每个键位尝试断键
    best_sim_overall = 0.0
    best_gen_overall = []
    best_tpl_overall = ""
    best_bond_overall = (0, 0)

    for bi in bonds:
        tpls = bi.get("templates", [])
        atoms = tuple(bi["atoms"])
        case.n_templates_matched += len(tpls)

        if not tpls:
            continue

        tpl_name = tpls[0]["name"]
        try:
            br = execute_disconnection(product, atoms, tpl_name)
        except Exception:
            continue

        attempt = {
            "bond": atoms,
            "template": tpl_name,
            "success": br.success,
            "is_fallback": br.is_fallback,
            "precursors": br.precursors if br.success else [],
            "error": br.error,
        }

        if not br.success or not br.precursors:
            case.all_attempts.append(attempt)
            continue

        case.n_bonds_succeeded += 1
        if br.is_fallback:
            case.used_fallback = True

        gen_set = {canonical(s) for s in br.precursors if canonical(s)}

        # 计算这次断键的最佳 tanimoto
        best_sim_this = 0.0
        for gen_smi in br.precursors:
            gen_mol = parse_mol(gen_smi)
            if gen_mol is None:
                continue
            for tm in true_mols:
                try:
                    sim = tanimoto(gen_mol, tm)
                    best_sim_this = max(best_sim_this, sim)
                except Exception:
                    pass

        attempt["best_tanimoto"] = round(best_sim_this, 4)
        attempt["exact_set"] = gen_set == true_set
        case.all_attempts.append(attempt)

        if best_sim_this > best_sim_overall:
            best_sim_overall = best_sim_this
            best_gen_overall = br.precursors
            best_tpl_overall = tpl_name
            best_bond_overall = atoms

        if gen_set == true_set:
            case.exact_set_match = True

    case.best_tanimoto = round(best_sim_overall, 4)
    case.best_gen_precursors = best_gen_overall
    case.best_template = best_tpl_overall
    case.best_bond = best_bond_overall

    # ===== FGI 模板尝试 =====
    try:
        scan_all = scan_applicable_reactions(product, mode="retro")
    except Exception:
        scan_all = {"ok": False}

    if scan_all.get("ok"):
        for tm in scan_all.get("matches", []):
            rxn = tm.rxn_smarts
            if not rxn or ">>" not in rxn:
                continue
            prod_side = rxn.split(">>")[-1]
            if "." in prod_side:
                continue  # 断键模板，跳过

            case.fgi_attempted += 1
            try:
                fgi_result = execute_fgi(product, tm.template_id)
            except Exception:
                continue

            if not fgi_result.success or not fgi_result.precursors:
                continue

            case.fgi_succeeded += 1

            fgi_sim = 0.0
            for gen_smi in fgi_result.precursors:
                gen_mol = parse_mol(gen_smi)
                if gen_mol is None:
                    continue
                for tm_mol in true_mols:
                    try:
                        sim = tanimoto(gen_mol, tm_mol)
                        fgi_sim = max(fgi_sim, sim)
                    except Exception:
                        pass

            if fgi_sim > case.best_fgi_tanimoto:
                case.best_fgi_tanimoto = round(fgi_sim, 4)
                case.best_fgi_precursors = fgi_result.precursors
                case.best_fgi_template = tm.template_id

    # 综合最佳
    case.combined_best_tanimoto = round(max(case.best_tanimoto, case.best_fgi_tanimoto), 4)

    # 失败分类 (基于综合最佳)
    best = case.combined_best_tanimoto
    if case.n_bonds_succeeded == 0 and case.fgi_succeeded == 0:
        case.failure_category = "all_failed"
    elif best < 0.3:
        case.failure_category = "very_low_similarity"
    elif best < 0.5:
        case.failure_category = "low_similarity"
    elif best < 0.7:
        case.failure_category = "moderate_similarity"
    elif best < 0.9:
        case.failure_category = "high_similarity_not_exact"
    elif not case.exact_set_match:
        case.failure_category = "near_miss"
    else:
        case.failure_category = "exact_match"

    return case


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------

def print_case_detail(case: CaseResult, rank: int):
    """打印单条案例的详细信息"""
    print(f"\n  {'─'*56}")
    print(f"  #{rank}  row={case.row_index}  tanimoto={case.combined_best_tanimoto}  [{case.failure_category}]")
    print(f"  {'─'*56}")
    print(f"  产物: {case.product[:90]}")
    print(f"  特征: {case.formula}  MW={case.mw}  rings={case.n_rings}  atoms={case.n_heavy}  hetero=[{case.heteroatoms}]")
    print(f"  真实前体: {case.true_precursors[:90]}")

    if case.best_gen_precursors:
        print(f"  最佳断键: {'.'.join(case.best_gen_precursors)[:90]}  (t={case.best_tanimoto})")
        print(f"    模板: {case.best_template}  键位: {case.best_bond}")
    else:
        print(f"  最佳断键: (无)")

    if case.best_fgi_precursors:
        print(f"  最佳FGI:  {'.'.join(case.best_fgi_precursors)[:90]}  (t={case.best_fgi_tanimoto})")
        print(f"    模板: {case.best_fgi_template}")

    print(f"  断键: 发现={case.n_bonds_found}  成功={case.n_bonds_succeeded}  "
          f"模板数={case.n_templates_matched}  fallback={'是' if case.used_fallback else '否'}")
    print(f"  FGI:  尝试={case.fgi_attempted}  成功={case.fgi_succeeded}")

    # 显示各断键尝试的 tanimoto
    successful = [a for a in case.all_attempts if a["success"]]
    if successful:
        print(f"  各键位 tanimoto:")
        for a in sorted(successful, key=lambda x: -x.get("best_tanimoto", 0))[:5]:
            fb = " [fallback]" if a.get("is_fallback") else ""
            prec_str = ".".join(a["precursors"])[:60]
            print(f"    bond={a['bond']}  t={a.get('best_tanimoto', 0):.4f}{fb}  → {prec_str}")


def print_summary(cases: List[CaseResult]):
    """打印汇总统计"""
    total = len(cases)
    if total == 0:
        print("  无数据")
        return

    # 按 failure_category 统计
    cat_counts = Counter(c.failure_category for c in cases)
    # tanimoto 分布
    tanimotos = [c.best_tanimoto for c in cases]
    tanimotos.sort()
    n = len(tanimotos)

    print(f"\n{'='*60}")
    print(f"  低 Tanimoto 案例分析汇总")
    print(f"{'='*60}")
    print(f"  总记录: {total}")
    print(f"  Tanimoto: mean={sum(tanimotos)/n:.4f}  median={tanimotos[n//2]:.4f}")
    print(f"            min={tanimotos[0]:.4f}  max={tanimotos[-1]:.4f}")

    print(f"\n  失败分类:")
    for cat, cnt in sorted(cat_counts.items(), key=lambda x: -x[1]):
        print(f"    {cat:30s}  {cnt:4d}  ({cnt/total:.1%})")

    # FGI 提升统计
    fgi_improved = [c for c in cases if c.best_fgi_tanimoto > c.best_tanimoto]
    fgi_any = [c for c in cases if c.fgi_succeeded > 0]
    print(f"\n  FGI 模板效果:")
    print(f"    有 FGI 匹配: {len(fgi_any)}/{total} ({len(fgi_any)/total:.1%})")
    print(f"    FGI 提升了 tanimoto: {len(fgi_improved)}/{total} ({len(fgi_improved)/total:.1%})")
    if fgi_improved:
        avg_boost = sum(c.best_fgi_tanimoto - c.best_tanimoto for c in fgi_improved) / len(fgi_improved)
        print(f"    平均提升幅度: +{avg_boost:.4f}")

    # 综合 vs 仅断键对比
    bond_only_tans = [c.best_tanimoto for c in cases]
    combined_tans = [c.combined_best_tanimoto for c in cases]
    print(f"\n  综合对比 (断键 vs 断键+FGI):")
    print(f"    仅断键 mean={sum(bond_only_tans)/n:.4f}  median={sorted(bond_only_tans)[n//2]:.4f}")
    print(f"    综合   mean={sum(combined_tans)/n:.4f}  median={sorted(combined_tans)[n//2]:.4f}")

    # 低分案例的产物特征统计
    low_cases = [c for c in cases if c.best_tanimoto < 0.5]
    if low_cases:
        print(f"\n  低分案例 (tanimoto < 0.5) 产物特征:")
        avg_rings = sum(c.n_rings for c in low_cases) / len(low_cases)
        avg_heavy = sum(c.n_heavy for c in low_cases) / len(low_cases)
        avg_mw = sum(c.mw for c in low_cases) / len(low_cases)
        print(f"    平均环数: {avg_rings:.1f}  平均重原子: {avg_heavy:.1f}  平均MW: {avg_mw:.0f}")

        # 杂原子分布
        hetero_counter = Counter()
        for c in low_cases:
            for h in c.heteroatoms.split(","):
                if h:
                    hetero_counter[h] += 1
        if hetero_counter:
            print(f"    杂原子分布: {dict(hetero_counter.most_common(10))}")

        # fallback 比例
        fb_count = sum(1 for c in low_cases if c.used_fallback)
        no_break = sum(1 for c in low_cases if c.n_bonds_succeeded == 0)
        print(f"    使用 fallback: {fb_count}/{len(low_cases)} ({fb_count/len(low_cases):.1%})")
        print(f"    所有断键失败: {no_break}/{len(low_cases)} ({no_break/len(low_cases):.1%})")

    # 高分但非精确匹配
    near = [c for c in cases if c.best_tanimoto >= 0.9 and not c.exact_set_match]
    if near:
        print(f"\n  Near-miss (tanimoto >= 0.9 但非精确匹配): {len(near)} 条")

    print(f"\n{'='*60}")


# ---------------------------------------------------------------------------
# HTML 输出
# ---------------------------------------------------------------------------

def generate_html(cases: List[CaseResult], file_label: str) -> str:
    """生成 HTML 可视化报告"""
    from rdkit.Chem.Draw import rdMolDraw2D

    def mol_svg(smiles: str, w: int = 300, h: int = 180) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return f'<span style="color:red">parse failed: {smiles[:40]}</span>'
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    parts = [f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>低Tanimoto案例分析 - {file_label}</title>
<style>
body {{ font-family: -apple-system, sans-serif; margin: 20px; background: #f5f5f5; }}
.card {{ background: white; border-radius: 8px; padding: 16px; margin: 12px 0;
         box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
.row {{ display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }}
.tag {{ padding: 2px 8px; border-radius: 4px; font-size: 12px; display: inline-block; margin: 2px; }}
.bad {{ background: #ffebee; color: #c62828; }}
.ok {{ background: #e8f5e9; color: #2e7d32; }}
.mid {{ background: #fff3e0; color: #e65100; }}
table {{ border-collapse: collapse; margin: 8px 0; font-size: 13px; }}
td, th {{ padding: 4px 10px; border: 1px solid #ddd; }}
th {{ background: #f5f5f5; text-align: left; }}
svg {{ max-width: 100%; }}
h2 {{ color: #333; border-bottom: 2px solid #1976d2; padding-bottom: 4px; }}
.smi {{ font-family: monospace; font-size: 11px; color: #666; word-break: break-all; }}
</style></head><body>
<h2>低 Tanimoto 案例分析 — {file_label} ({len(cases)} 条)</h2>
"""]

    for i, c in enumerate(cases):
        t_cls = "bad" if c.combined_best_tanimoto < 0.5 else ("mid" if c.combined_best_tanimoto < 0.9 else "ok")
        parts.append(f'<div class="card">')
        parts.append(f'<b>#{i+1}</b> row={c.row_index} '
                     f'<span class="tag {t_cls}">combined={c.combined_best_tanimoto}</span> '
                     f'<span class="tag">断键={c.best_tanimoto} FGI={c.best_fgi_tanimoto}</span> '
                     f'<span class="tag">{c.failure_category}</span> '
                     f'<span class="tag">{c.formula} MW={c.mw}</span>')

        parts.append('<div class="row" style="margin-top:8px">')
        # 产物
        parts.append(f'<div style="text-align:center"><div style="font-size:11px;color:#888">产物</div>{mol_svg(c.product)}'
                     f'<div class="smi">{c.product[:80]}</div></div>')

        # 真实前体
        for j, frag in enumerate(c.true_precursors.split(".")):
            if j == 0:
                parts.append('<div style="font-size:24px;color:#666;padding:0 4px">⇐</div>')
            else:
                parts.append('<div style="font-size:18px;color:#999">+</div>')
            parts.append(f'<div style="text-align:center"><div style="font-size:11px;color:#888">真实前体{j+1}</div>'
                         f'{mol_svg(frag, 250, 150)}<div class="smi">{frag[:60]}</div></div>')
        parts.append('</div>')

        # 最佳生成前体
        if c.best_gen_precursors:
            parts.append(f'<div style="margin-top:8px;font-size:12px;color:#555">'
                         f'最佳断键 (t={c.best_tanimoto}, 模板: {c.best_template}, 键位: {c.best_bond}):</div>')
            parts.append('<div class="row">')
            for j, gen in enumerate(c.best_gen_precursors):
                if j > 0:
                    parts.append('<div style="font-size:18px;color:#999">+</div>')
                parts.append(f'<div style="text-align:center">{mol_svg(gen, 250, 150)}'
                             f'<div class="smi">{gen[:60]}</div></div>')
            parts.append('</div>')

        # FGI 最佳前体
        if c.best_fgi_precursors:
            fgi_cls = "ok" if c.best_fgi_tanimoto > c.best_tanimoto else "mid"
            parts.append(f'<div style="margin-top:8px;font-size:12px;color:#1976d2">'
                         f'最佳FGI (t={c.best_fgi_tanimoto}, 模板: {c.best_fgi_template}):</div>')
            parts.append('<div class="row">')
            for j, gen in enumerate(c.best_fgi_precursors):
                if j > 0:
                    parts.append('<div style="font-size:18px;color:#999">+</div>')
                parts.append(f'<div style="text-align:center">{mol_svg(gen, 250, 150)}'
                             f'<div class="smi">{gen[:60]}</div></div>')
            parts.append('</div>')

        parts.append('</div>')

    parts.append('</body></html>')
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="低 Tanimoto 案例分析")
    parser.add_argument("--file", default="BTF1", help="数据文件 (BTF1-10, ATF, NEW)")
    parser.add_argument("--size", type=int, default=200, help="采样数量")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    parser.add_argument("--top", type=int, default=20, help="显示最差的 N 条")
    parser.add_argument("--threshold", type=float, default=None,
                        help="只显示 tanimoto < threshold 的案例")
    parser.add_argument("--html", action="store_true", help="输出 HTML 可视化")
    args = parser.parse_args()

    file_key = args.file.upper()
    if file_key not in FILE_MAP:
        print(f"未知文件: {args.file}，可选: {', '.join(FILE_MAP.keys())}")
        sys.exit(1)

    path = FILE_MAP[file_key]
    if not path.exists():
        print(f"文件不存在: {path}")
        sys.exit(1)

    records = load_csv(str(path), sample_size=args.size, random_seed=args.seed)
    print(f"加载 {file_key}: {len(records)} 条")

    # 评估所有记录
    print("评估中...", flush=True)
    t0 = time.time()
    cases = []
    for i, rec in enumerate(records):
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(records)}...", flush=True)
        cases.append(evaluate_single(rec))
    elapsed = time.time() - t0
    print(f"评估完成: {elapsed:.1f}s")

    # 汇总
    print_summary(cases)

    # 筛选要展示的案例
    if args.threshold is not None:
        show_cases = [c for c in cases if c.combined_best_tanimoto < args.threshold]
        show_cases.sort(key=lambda c: c.combined_best_tanimoto)
        label = f"tanimoto < {args.threshold}"
    else:
        show_cases = [c for c in cases if not c.exact_set_match]
        show_cases.sort(key=lambda c: c.combined_best_tanimoto)
        show_cases = show_cases[:args.top]
        label = f"最差 {args.top} 条"

    print(f"\n  展示: {label} ({len(show_cases)} 条)")
    for i, c in enumerate(show_cases):
        print_case_detail(c, i + 1)

    # HTML 输出
    if args.html:
        html = generate_html(show_cases, f"{file_key} - {label}")
        out_path = Path(f"near_miss_{file_key.lower()}.html")
        out_path.write_text(html, encoding="utf-8")
        print(f"\nHTML 输出: {out_path.resolve()}")
        import webbrowser
        webbrowser.open(str(out_path.resolve()))


if __name__ == "__main__":
    main()
