#!/usr/bin/env python
"""LLM 逆合成决策平台 — 为 LLM 提供完整的分子分析和断键候选方案。

基于 llm_bond_test.py 改造，从交互式工具升级为：
  1. 单分子模式: 输入 SMILES，输出完整的 LLM 决策上下文
  2. 批量验证模式: 从数据集读取，对比 LLM 选择 vs 真实前体

用法:
    # 单分子分析 — 输出 LLM 可读的决策上下文
    python Rachel/tools/llm_retro_platform.py --smiles "c1ccc(-c2ccccc2)cc1"

    # 批量验证 — 从数据集读取，输出每条记录的决策上下文
    python Rachel/tools/llm_retro_platform.py --file BTF1 --size 10

    # JSON 输出 — 方便程序化处理
    python Rachel/tools/llm_retro_platform.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --json

输出内容:
  - 分子基本信息 (结构、官能团、复杂度)
  - 所有可断键位 + 每个键位的所有等价断法预览
  - FGI 模板匹配
  - 官能团兼容性警告
  - (批量模式) 真实前体对比
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.chem_tools._rdkit_utils import canonical, parse_mol, tanimoto
from Rachel.chem_tools.mol_info import analyze_molecule
from Rachel.chem_tools.fg_detect import detect_functional_groups
from Rachel.chem_tools.cs_score import compute_cs_score, classify_complexity
from Rachel.chem_tools.template_scan import (
    find_disconnectable_bonds,
    scan_applicable_reactions,
)
from Rachel.chem_tools.bond_break import (
    execute_disconnection,
    preview_disconnections,
    execute_fgi,
)
from Rachel.chem_tools.fg_warnings import check_fg_conflicts
from Rachel.chem_tools.smart_cap import suggest_capping


DATA_DIR = ROOT / "Rachel" / "data"
FILE_MAP = {
    "ATF":  DATA_DIR / "USPTO50K" / "datasetATF.csv",
    **{f"BTF{i}": DATA_DIR / "USPTO50K" / f"datasetBTF{i}.csv" for i in range(1, 11)},
    "NEW":  DATA_DIR / "NEW" / "usptonew.csv",
}


# ─────────────────────────────────────────────────────────────────────────
# 核心: 生成 LLM 决策上下文
# ─────────────────────────────────────────────────────────────────────────

def build_decision_context(smiles: str) -> Dict[str, Any]:
    """为一个目标分子生成完整的 LLM 逆合成决策上下文。

    Returns dict with keys:
      - molecule: 基本分子信息
      - functional_groups: 官能团列表
      - complexity: CS 评分和分类
      - disconnectable_bonds: 可断键位 + 每个键位的所有等价方案预览
      - fgi_options: FGI 模板匹配
      - warnings: 官能团兼容性警告
    """
    ctx: Dict[str, Any] = {"smiles": canonical(smiles) or smiles}

    # ── 1. 分子信息 ──
    mol_info = analyze_molecule(smiles)
    descriptors = mol_info.get("descriptors", {})
    stereo_info = mol_info.get("stereo", {})
    ctx["molecule"] = {
        "canonical": mol_info.get("smiles", ""),
        "formula": mol_info.get("formula", ""),
        "mw": mol_info.get("mw", 0),
        "heavy_atoms": len(mol_info.get("atoms", [])),
        "rings": descriptors.get("ring_count", 0),
        "aromatic_rings": descriptors.get("aromatic_ring_count", 0),
        "stereocenters": len(stereo_info.get("chiral_centers", [])),
        "rotatable_bonds": descriptors.get("rotatable_bonds", 0),
    }

    # ── 2. 官能团 ──
    fg = detect_functional_groups(smiles)
    fg_list = []
    groups_data = fg.get("groups", {})
    if isinstance(groups_data, dict):
        for fg_name, fg_info in groups_data.items():
            fg_list.append({
                "name": fg_name,
                "count": fg_info.get("count", 1) if isinstance(fg_info, dict) else 1,
                "atoms": fg_info.get("atoms", []) if isinstance(fg_info, dict) else [],
            })
    ctx["functional_groups"] = fg_list

    # ── 3. 复杂度 ──
    cs = compute_cs_score(smiles)
    cls = classify_complexity(smiles)
    ctx["complexity"] = {
        "cs_score": cs.get("cs_score", 0),
        "classification": cls.get("classification", ""),
        "is_terminal": cls.get("classification", "") == "trivial",
        "dimensions": cs.get("dimensions", {}),
    }

    # ── 4. 可断键位 + 等价方案预览 ──
    scan = find_disconnectable_bonds(smiles)
    bonds_ctx = []
    if scan.get("ok") and scan.get("bonds"):
        for bond_info in scan["bonds"]:
            atoms = tuple(bond_info["atoms"])
            bond_entry = {
                "atoms": list(atoms),
                "bond_type": bond_info.get("bond_type", ""),
                "in_ring": bond_info.get("in_ring", False),
                "heuristic_score": bond_info.get("heuristic_score", 0),
                "n_templates": len(bond_info.get("templates", [])),
                "template_names": [
                    t.get("name", "") for t in bond_info.get("templates", [])
                ],
            }

            # 预览所有等价断法
            preview = preview_disconnections(smiles, atoms)
            if preview.get("ok"):
                alts = []
                for alt in preview.get("alternatives", []):
                    alts.append({
                        "template": alt["name"],
                        "template_id": alt["template_id"],
                        "precursors": alt["precursors"],
                        "incompatible_with": alt.get("incompatible_groups", []),
                    })
                bond_entry["alternatives"] = alts
            else:
                bond_entry["alternatives"] = []

            # Smart capping 推理（模板无关的规则推断）
            try:
                cap_result = suggest_capping(smiles, tuple(atoms))
                if cap_result.get("ok") and cap_result.get("proposals"):
                    bond_entry["smart_capping"] = cap_result["proposals"][:3]
            except Exception:
                pass

            bonds_ctx.append(bond_entry)

    ctx["disconnectable_bonds"] = bonds_ctx
    ctx["n_bonds"] = len(bonds_ctx)

    # ── 5. FGI 选项 ──
    fgi_ctx = []
    try:
        fgi_scan = scan_applicable_reactions(smiles, mode="retro")
        if fgi_scan.get("ok"):
            for tm in fgi_scan.get("matches", []):
                rxn = tm.rxn_smarts
                if not rxn or ">>" not in rxn:
                    continue
                prod_side = rxn.split(">>")[-1]
                if "." in prod_side:
                    continue  # 断键模板，跳过
                try:
                    fgi_r = execute_fgi(smiles, tm.template_id)
                    if fgi_r.success and fgi_r.precursors:
                        fgi_ctx.append({
                            "template": tm.name,
                            "template_id": tm.template_id,
                            "precursors": fgi_r.precursors,
                        })
                except Exception:
                    pass
    except Exception:
        pass
    ctx["fgi_options"] = fgi_ctx

    # ── 6. 官能团警告 ──
    warnings = check_fg_conflicts(smiles)
    ctx["warnings"] = warnings.get("conflicts", []) if isinstance(warnings, dict) else []

    return ctx


# ─────────────────────────────────────────────────────────────────────────
# 格式化输出: 人类可读 / LLM prompt 友好
# ─────────────────────────────────────────────────────────────────────────

def format_context_text(ctx: Dict[str, Any], true_precursors: str = "") -> str:
    """将决策上下文格式化为 LLM 可读的文本。"""
    lines = []
    mol = ctx.get("molecule", {})
    lines.append(f"═══ 目标分子 ═══")
    lines.append(f"SMILES:    {ctx['smiles']}")
    lines.append(f"分子式:    {mol.get('formula', '?')}  MW={mol.get('mw', 0)}")
    lines.append(f"重原子:    {mol.get('heavy_atoms', 0)}  环={mol.get('rings', 0)}  "
                 f"手性中心={mol.get('stereocenters', 0)}")

    # 复杂度
    cplx = ctx.get("complexity", {})
    lines.append(f"复杂度:    CS={cplx.get('cs_score', 0):.2f}  "
                 f"分类={cplx.get('classification', '?')}  "
                 f"{'(可作为终端)' if cplx.get('is_terminal') else ''}")

    # 官能团
    fgs = ctx.get("functional_groups", [])
    if fgs:
        fg_names = [f"{g['name']}(×{g['count']})" if g['count'] > 1 else g['name']
                    for g in fgs[:15]]
        lines.append(f"官能团:    {', '.join(fg_names)}")

    # 警告
    warns = ctx.get("warnings", [])
    if warns:
        lines.append(f"⚠ 兼容性警告: {len(warns)} 条")
        for w in warns[:5]:
            if isinstance(w, dict):
                lines.append(f"  - {w.get('message', w)}")
            else:
                lines.append(f"  - {w}")

    # 真实前体 (验证模式)
    if true_precursors:
        lines.append(f"\n── 真实前体 (参考答案) ──")
        lines.append(f"  {true_precursors}")

    # 可断键位
    bonds = ctx.get("disconnectable_bonds", [])
    lines.append(f"\n═══ 可断键位: {len(bonds)} 个 ═══")
    for i, b in enumerate(bonds):
        alts = b.get("alternatives", [])
        lines.append(f"\n  键位 #{i}: atoms={b['atoms']}  type={b['bond_type']}  "
                     f"ring={'是' if b['in_ring'] else '否'}  "
                     f"score={b['heuristic_score']:.4f}  "
                     f"方案数={len(alts)}")

        if alts:
            for j, alt in enumerate(alts):
                prec_str = " + ".join(alt["precursors"])
                incompat = ""
                if alt.get("incompatible_with"):
                    incompat = f"  [不兼容: {', '.join(alt['incompatible_with'][:3])}]"
                lines.append(f"    [{j}] {alt['template']}")
                lines.append(f"        → {prec_str}{incompat}")
        else:
            for tname in b.get("template_names", []):
                lines.append(f"    - {tname}")

    # FGI 选项
    fgi = ctx.get("fgi_options", [])
    if fgi:
        lines.append(f"\n═══ FGI 选项: {len(fgi)} 个 ═══")
        for i, f in enumerate(fgi[:10]):
            prec_str = " + ".join(f["precursors"])
            lines.append(f"  [{i}] {f['template']}")
            lines.append(f"      → {prec_str}")

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────
# 批量验证: 对比工具层 vs 真实前体
# ─────────────────────────────────────────────────────────────────────────

def evaluate_with_context(
    smiles: str,
    true_precursors: str,
    ctx: Dict[str, Any],
) -> Dict[str, Any]:
    """用决策上下文评估: 工具层的候选方案中是否包含真实前体。"""
    true_frags = true_precursors.split(".")
    true_set = {canonical(f) for f in true_frags if canonical(f)}
    true_mols = [m for f in true_frags if (m := parse_mol(f))]

    result = {
        "product": smiles,
        "true_precursors": true_precursors,
        "n_bonds": ctx.get("n_bonds", 0),
        "best_tanimoto": 0.0,
        "exact_match": False,
        "best_template": "",
        "best_precursors": [],
    }

    if not true_set or not true_mols:
        return result

    # 检查所有断键方案
    for bond in ctx.get("disconnectable_bonds", []):
        for alt in bond.get("alternatives", []):
            gen_set = {canonical(s) for s in alt["precursors"] if canonical(s)}
            if gen_set == true_set:
                result["exact_match"] = True
                result["best_tanimoto"] = 1.0
                result["best_template"] = alt["template"]
                result["best_precursors"] = alt["precursors"]
                return result

            # Tanimoto
            for gen_smi in alt["precursors"]:
                gen_mol = parse_mol(gen_smi)
                if gen_mol is None:
                    continue
                for tm in true_mols:
                    try:
                        sim = tanimoto(gen_mol, tm)
                        if sim > result["best_tanimoto"]:
                            result["best_tanimoto"] = sim
                            result["best_template"] = alt["template"]
                            result["best_precursors"] = alt["precursors"]
                    except Exception:
                        pass

    # 检查 FGI 方案
    for fgi in ctx.get("fgi_options", []):
        gen_set = {canonical(s) for s in fgi["precursors"] if canonical(s)}
        if gen_set == true_set:
            result["exact_match"] = True
            result["best_tanimoto"] = 1.0
            result["best_template"] = fgi["template"]
            result["best_precursors"] = fgi["precursors"]
            return result

        for gen_smi in fgi["precursors"]:
            gen_mol = parse_mol(gen_smi)
            if gen_mol is None:
                continue
            for tm in true_mols:
                try:
                    sim = tanimoto(gen_mol, tm)
                    if sim > result["best_tanimoto"]:
                        result["best_tanimoto"] = sim
                        result["best_template"] = fgi["template"]
                        result["best_precursors"] = fgi["precursors"]
                except Exception:
                    pass

    return result


# ─────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="LLM 逆合成决策平台")
    parser.add_argument("--smiles", type=str, default=None,
                        help="单分子 SMILES")
    parser.add_argument("--file", type=str, default=None,
                        help="数据文件 (BTF1-10, ATF, NEW)")
    parser.add_argument("--size", type=int, default=10,
                        help="批量模式采样数量")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子")
    parser.add_argument("--json", action="store_true",
                        help="输出 JSON 格式")
    parser.add_argument("--index", type=int, default=None,
                        help="批量模式中只处理第 N 条记录 (0-based)")
    args = parser.parse_args()

    if args.smiles:
        # ── 单分子模式 ──
        mol = parse_mol(args.smiles)
        if mol is None:
            print(f"无效 SMILES: {args.smiles}")
            sys.exit(1)

        t0 = time.time()
        ctx = build_decision_context(args.smiles)
        elapsed = time.time() - t0

        if args.json:
            ctx["_elapsed_sec"] = round(elapsed, 2)
            print(json.dumps(ctx, indent=2, ensure_ascii=False))
        else:
            print(format_context_text(ctx))
            print(f"\n  (耗时: {elapsed:.1f}s)")

    elif args.file:
        # ── 批量验证模式 ──
        from Rachel.tests.data_driven.loader import load_csv

        file_key = args.file.upper()
        if file_key not in FILE_MAP:
            print(f"未知文件: {args.file}，可选: {', '.join(FILE_MAP.keys())}")
            sys.exit(1)

        path = FILE_MAP[file_key]
        if not path.exists():
            print(f"文件不存在: {path}")
            sys.exit(1)

        records = load_csv(str(path), sample_size=args.size, random_seed=args.seed)
        print(f"加载 {file_key}: {len(records)} 条\n")

        if args.index is not None:
            if args.index < 0 or args.index >= len(records):
                print(f"索引超出范围: {args.index} (共 {len(records)} 条)")
                sys.exit(1)
            records = [records[args.index]]

        exact_count = 0
        total = len(records)
        all_results = []

        for i, rec in enumerate(records):
            t0 = time.time()
            ctx = build_decision_context(rec.product_smiles)
            elapsed = time.time() - t0

            ev = evaluate_with_context(
                rec.product_smiles, rec.precursors_smiles, ctx
            )
            ev["_elapsed_sec"] = round(elapsed, 2)
            ev["_row_index"] = rec.row_index
            all_results.append(ev)

            if ev["exact_match"]:
                exact_count += 1

            if args.json:
                continue

            # 文本输出
            status = "✓ EXACT" if ev["exact_match"] else f"  t={ev['best_tanimoto']:.4f}"
            print(f"[{i+1}/{total}] row={rec.row_index}  {status}  "
                  f"bonds={ev['n_bonds']}  ({elapsed:.1f}s)")

            if not ev["exact_match"]:
                print(f"  产物:   {rec.product_smiles[:80]}")
                print(f"  真实:   {rec.precursors_smiles[:80]}")
                if ev["best_precursors"]:
                    print(f"  最佳:   {'.'.join(ev['best_precursors'])[:80]}")
                    print(f"  模板:   {ev['best_template']}")

            # 如果是单条模式 (--index)，输出完整上下文
            if args.index is not None:
                print()
                print(format_context_text(ctx, rec.precursors_smiles))

        if args.json:
            print(json.dumps(all_results, indent=2, ensure_ascii=False))
        else:
            print(f"\n{'='*50}")
            print(f"  exact_match: {exact_count}/{total} = {exact_count/total:.1%}")
            print(f"{'='*50}")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
