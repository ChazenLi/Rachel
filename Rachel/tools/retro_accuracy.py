"""严格逆合成准确率评估 — 三个级别的匹配标准。

用法:
    python Rachel/tools/retro_accuracy.py                    # BTF1, 200条
    python Rachel/tools/retro_accuracy.py --file BTF1 --size 500
    python Rachel/tools/retro_accuracy.py --file ALL --size 100

三个匹配级别:
  - exact_set:  生成前体集合 == 真实前体集合 (canonical SMILES 集合完全一致)
  - partial:    生成前体集合与真实前体集合有交集 (至少一个片段匹配)
  - tanimoto:   最佳 Tanimoto 相似度统计
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

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


def canonicalize_set(smiles_list: List[str]) -> Set[str]:
    """将 SMILES 列表转为 canonical 集合"""
    result = set()
    for s in smiles_list:
        c = canonical(s)
        if c:
            result.add(c)
    return result


def evaluate_records(records, verbose: bool = False):
    """对每条记录做严格的逆合成评估（断键 + FGI）"""

    total = 0
    no_bonds = 0
    all_breaks_failed = 0

    # ── 仅断键 ──
    exact_set_matches = 0
    partial_matches = 0
    best_tanimotos: List[float] = []
    bonds_tried_total = 0
    bonds_succeeded_total = 0

    # ── 断键 + FGI 综合 ──
    combined_exact = 0
    combined_partial = 0
    combined_tanimotos: List[float] = []
    fgi_attempted_total = 0
    fgi_succeeded_total = 0
    fgi_improved_count = 0  # FGI 比断键更好的记录数

    for rec in records:
        total += 1
        product = rec.product_smiles

        # 真实前体
        true_frags = rec.precursors_smiles.split(".")
        true_set = canonicalize_set(true_frags)
        true_mols = [m for f in true_frags if (m := parse_mol(f))]

        if not true_set:
            continue

        # ===== 断键路径 =====
        best_sim_bond = 0.0
        found_exact_bond = False
        found_partial_bond = False
        any_bond_success = False

        try:
            scan = find_disconnectable_bonds(product)
        except Exception:
            scan = {"ok": False}

        if scan.get("ok") and scan.get("bonds"):
            for bi in scan["bonds"]:
                tpls = bi.get("templates", [])
                if not tpls:
                    continue
                bonds_tried_total += 1

                try:
                    br = execute_disconnection(product, tuple(bi["atoms"]), tpls[0]["name"])
                except Exception:
                    continue

                if not br.success or not br.precursors:
                    continue

                bonds_succeeded_total += 1
                any_bond_success = True

                gen_set = canonicalize_set(br.precursors)
                if gen_set == true_set:
                    found_exact_bond = True
                if gen_set & true_set:
                    found_partial_bond = True

                for gen_smi in br.precursors:
                    gen_mol = parse_mol(gen_smi)
                    if gen_mol is None:
                        continue
                    for tm in true_mols:
                        try:
                            sim = tanimoto(gen_mol, tm)
                            best_sim_bond = max(best_sim_bond, sim)
                        except Exception:
                            pass
        else:
            no_bonds += 1

        # ===== FGI 路径 =====
        best_sim_fgi = 0.0
        found_exact_fgi = False
        found_partial_fgi = False
        any_fgi_success = False

        try:
            scan_all = scan_applicable_reactions(product, mode="retro")
        except Exception:
            scan_all = {"ok": False}

        if scan_all.get("ok"):
            for tm_match in scan_all.get("matches", []):
                rxn = tm_match.rxn_smarts
                if not rxn or ">>" not in rxn:
                    continue
                prod_side = rxn.split(">>")[-1]
                if "." in prod_side:
                    continue  # 断键模板，跳过

                fgi_attempted_total += 1
                try:
                    fgi_r = execute_fgi(product, tm_match.template_id)
                except Exception:
                    continue

                if not fgi_r.success or not fgi_r.precursors:
                    continue

                fgi_succeeded_total += 1
                any_fgi_success = True

                gen_set = canonicalize_set(fgi_r.precursors)
                if gen_set == true_set:
                    found_exact_fgi = True
                if gen_set & true_set:
                    found_partial_fgi = True

                for gen_smi in fgi_r.precursors:
                    gen_mol = parse_mol(gen_smi)
                    if gen_mol is None:
                        continue
                    for tm in true_mols:
                        try:
                            sim = tanimoto(gen_mol, tm)
                            best_sim_fgi = max(best_sim_fgi, sim)
                        except Exception:
                            pass

        # ===== 汇总 =====
        any_success = any_bond_success or any_fgi_success

        if not any_bond_success and not any_fgi_success:
            if not scan.get("ok") or not scan.get("bonds"):
                pass  # already counted in no_bonds
            else:
                all_breaks_failed += 1
            continue

        # 仅断键统计
        if any_bond_success:
            best_tanimotos.append(best_sim_bond)
            if found_exact_bond:
                exact_set_matches += 1
            if found_partial_bond:
                partial_matches += 1

        # 综合统计
        combined_sim = max(best_sim_bond, best_sim_fgi)
        combined_tanimotos.append(combined_sim)
        if found_exact_bond or found_exact_fgi:
            combined_exact += 1
        if found_partial_bond or found_partial_fgi:
            combined_partial += 1
        if best_sim_fgi > best_sim_bond:
            fgi_improved_count += 1

    # ── 输出 ──
    bond_evaluated = len(best_tanimotos)
    combined_evaluated = len(combined_tanimotos)
    best_tanimotos.sort()
    combined_tanimotos.sort()

    print(f"\n{'='*60}")
    print(f"  严格逆合成准确率评估 (断键 + FGI)")
    print(f"{'='*60}")
    print(f"  总记录数:           {total}")
    print(f"  无可断键位:         {no_bonds}")
    print(f"  所有方法失败:       {all_breaks_failed}")
    print(f"  断键成功评估:       {bond_evaluated}")
    print(f"  综合成功评估:       {combined_evaluated}")
    print(f"  尝试键位总数:       {bonds_tried_total}")
    print(f"  成功断键总数:       {bonds_succeeded_total}")
    print(f"  FGI 尝试/成功:      {fgi_attempted_total}/{fgi_succeeded_total}")
    print(f"  FGI 提升记录数:     {fgi_improved_count}")

    def _print_block(label, evaluated, exact, partial, tans):
        print(f"\n  ┌─ {label} (基于 {evaluated} 条) ─┐")
        print(f"  │")
        if evaluated:
            print(f"  │  exact_set: {exact}/{evaluated} = {exact/evaluated:.1%}")
            print(f"  │  partial:   {partial}/{evaluated} = {partial/evaluated:.1%}")
        else:
            print(f"  │  (无数据)")
        if tans:
            n = len(tans)
            mean_t = sum(tans) / n
            print(f"  │  Tanimoto:  mean={mean_t:.4f}  median={tans[n//2]:.4f}")
            print(f"  │             P25={tans[n//4]:.4f}   P75={tans[3*n//4]:.4f}")
            labels = ["<0.3", "0.3-0.5", "0.5-0.7", "0.7-0.9", "0.9-1.0", "=1.0"]
            counts = [0] * 6
            for t in tans:
                if t >= 1.0:    counts[5] += 1
                elif t >= 0.9:  counts[4] += 1
                elif t >= 0.7:  counts[3] += 1
                elif t >= 0.5:  counts[2] += 1
                elif t >= 0.3:  counts[1] += 1
                else:           counts[0] += 1
            print(f"  │  分布:")
            for lb, cnt in zip(labels, counts):
                bar = "█" * int(cnt / n * 30)
                print(f"  │    {lb:8s} {cnt:4d} ({cnt/n:5.1%}) {bar}")
        print(f"  │")
        print(f"  └{'─'*48}┘")

    _print_block("仅断键", bond_evaluated, exact_set_matches, partial_matches, best_tanimotos)
    _print_block("断键 + FGI 综合", combined_evaluated, combined_exact, combined_partial, combined_tanimotos)
    print()


def main():
    parser = argparse.ArgumentParser(description="严格逆合成准确率评估")
    parser.add_argument("--file", default="BTF1", help="数据文件 (BTF1-10, ATF, NEW, ALL)")
    parser.add_argument("--size", type=int, default=200, help="采样数量")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    args = parser.parse_args()

    file_key = args.file.upper()
    if file_key == "ALL":
        files = list(FILE_MAP.items())
    elif file_key in FILE_MAP:
        files = [(file_key, FILE_MAP[file_key])]
    else:
        print(f"未知文件: {args.file}")
        sys.exit(1)

    for label, path in files:
        if not path.exists():
            print(f"跳过 {label}: 文件不存在")
            continue
        records = load_csv(str(path), sample_size=args.size, random_seed=args.seed)
        print(f"\n加载 {label}: {len(records)} 条")
        t0 = time.time()
        evaluate_records(records)
        print(f"  耗时: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
