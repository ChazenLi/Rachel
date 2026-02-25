"""快速数据驱动验证脚本 — 独立于 pytest，直接跑化学层验证并输出报告。

用法:
    python Rachel/tools/run_validation.py                          # 默认 datasetBTF1, 100条
    python Rachel/tools/run_validation.py --file BTF3 --size 500   # 指定分片和数量
    python Rachel/tools/run_validation.py --file ALL --size 50     # 所有分片各跑50条
    python Rachel/tools/run_validation.py --file NEW --size 200    # 跑 usptonew.csv
    python Rachel/tools/run_validation.py --step retro             # 只跑逆合成评估
    python Rachel/tools/run_validation.py --step parse,bond        # 跑解析+断键
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# 确保项目根目录在 sys.path
ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# 静默 RDKit 警告（模板执行时大量 aromatic/kekulize 警告）— 必须在任何 chem_tools import 之前
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.chem_tools._rdkit_utils import canonical, parse_mol, tanimoto
from Rachel.chem_tools.bond_break import execute_disconnection
from Rachel.chem_tools.forward_validate import validate_forward
from Rachel.chem_tools.mol_info import analyze_molecule
from Rachel.chem_tools.template_scan import find_disconnectable_bonds
from Rachel.tests.data_driven.loader import load_csv, get_skipped_count


# ── 数据文件映射 ──────────────────────────────────────────────────
DATA_DIR = ROOT / "Rachel" / "data"
FILE_MAP = {
    "ATF":  DATA_DIR / "USPTO50K" / "datasetATF.csv",
    "BTF1": DATA_DIR / "USPTO50K" / "datasetBTF1.csv",
    "BTF2": DATA_DIR / "USPTO50K" / "datasetBTF2.csv",
    "BTF3": DATA_DIR / "USPTO50K" / "datasetBTF3.csv",
    "BTF4": DATA_DIR / "USPTO50K" / "datasetBTF4.csv",
    "BTF5": DATA_DIR / "USPTO50K" / "datasetBTF5.csv",
    "BTF6": DATA_DIR / "USPTO50K" / "datasetBTF6.csv",
    "BTF7": DATA_DIR / "USPTO50K" / "datasetBTF7.csv",
    "BTF8": DATA_DIR / "USPTO50K" / "datasetBTF8.csv",
    "BTF9": DATA_DIR / "USPTO50K" / "datasetBTF9.csv",
    "BTF10": DATA_DIR / "USPTO50K" / "datasetBTF10.csv",
    "NEW":  DATA_DIR / "NEW" / "usptonew.csv",
}

STEPS_ALL = ["parse", "bond", "forward", "pipeline", "retro"]


# ── 结果收集 ─────────────────────────────────────────────────────
@dataclass
class Stats:
    success: int = 0
    fail: int = 0
    errors: List[Dict[str, Any]] = field(default_factory=list)
    timings: List[float] = field(default_factory=list)

    @property
    def total(self) -> int:
        return self.success + self.fail

    @property
    def rate(self) -> float:
        return self.success / self.total if self.total else 0.0

    @property
    def avg_ms(self) -> float:
        return (sum(self.timings) / len(self.timings) * 1000) if self.timings else 0.0


def mol_features(smiles: str) -> Dict[str, Any]:
    mol = parse_mol(smiles)
    if not mol:
        return {}
    from rdkit.Chem import Descriptors
    return {
        "rings": mol.GetRingInfo().NumRings(),
        "hetero": sorted(set(a.GetSymbol() for a in mol.GetAtoms() if a.GetSymbol() != "C")),
        "mw": round(Descriptors.MolWt(mol)),
        "atoms": mol.GetNumHeavyAtoms(),
    }


# ── Step 1: SMILES 解析 ─────────────────────────────────────────
def run_parse(records, stats: Dict[str, Stats]):
    """验证 SMILES 可解析、规范化幂等、Round-Trip"""
    s = stats["parse"]
    for rec in records:
        ok = True
        # 产物可解析
        if parse_mol(rec.product_smiles) is None:
            ok = False
            s.errors.append({"smiles": rec.product_smiles, "reason": "product parse failed", "row": rec.row_index})
        # 前体各片段可解析
        for frag in rec.precursors_smiles.split("."):
            if parse_mol(frag) is None:
                ok = False
                s.errors.append({"smiles": frag, "reason": "precursor fragment parse failed", "row": rec.row_index})
        # 规范化幂等
        for smi in [rec.product_smiles] + rec.precursors_smiles.split("."):
            c1 = canonical(smi)
            if c1 and canonical(c1) != c1:
                ok = False
                s.errors.append({"smiles": smi, "reason": f"canonical not idempotent", "row": rec.row_index})
        if ok:
            s.success += 1
        else:
            s.fail += 1


# ── Step 2: 断键验证 ────────────────────────────────────────────
def run_bond(records, stats: Dict[str, Stats]):
    """验证键位发现 + 断键执行（前体可解析、重原子守恒）"""
    s = stats["bond"]
    for rec in records:
        product = rec.product_smiles
        product_mol = parse_mol(product)
        if product_mol is None:
            s.fail += 1
            continue

        t0 = time.time()
        try:
            scan = find_disconnectable_bonds(product)
        except Exception as e:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": f"scan exception: {e}", "row": rec.row_index})
            continue

        if not scan.get("ok") or not scan.get("bonds"):
            s.success += 1  # 没有可断键不算失败
            s.timings.append(time.time() - t0)
            continue

        product_heavy = product_mol.GetNumHeavyAtoms()
        rec_ok = True

        for bond_info in scan["bonds"]:
            templates = bond_info.get("templates", [])
            if not templates:
                continue
            atoms = bond_info["atoms"]
            tpl_name = templates[0]["name"]
            try:
                br = execute_disconnection(product, tuple(atoms), tpl_name)
            except Exception as e:
                rec_ok = False
                s.errors.append({"smiles": product, "reason": f"break exception: {e}", "bond": atoms, "row": rec.row_index})
                continue

            if not br.success:
                if not br.error:
                    rec_ok = False
                    s.errors.append({"smiles": product, "reason": "failed break has no error", "bond": atoms, "row": rec.row_index})
                continue

            # 前体可解析
            for p in br.precursors:
                if parse_mol(p) is None:
                    rec_ok = False
                    s.errors.append({"smiles": product, "reason": f"unparseable precursor: {p}", "bond": atoms, "row": rec.row_index})

            # 重原子守恒
            prec_heavy = sum(m.GetNumHeavyAtoms() for p in br.precursors if (m := parse_mol(p)))
            if prec_heavy < product_heavy:
                rec_ok = False
                s.errors.append({
                    "smiles": product, "reason": f"heavy atom loss: {product_heavy} -> {prec_heavy}",
                    "bond": atoms, "row": rec.row_index,
                })

        s.timings.append(time.time() - t0)
        if rec_ok:
            s.success += 1
        else:
            s.fail += 1


# ── Step 3: 正向验证 ────────────────────────────────────────────
def run_forward(records, stats: Dict[str, Stats]):
    """验证 validate_forward 输出结构、feasibility_score 范围"""
    s = stats["forward"]
    for rec in records:
        product = rec.product_smiles
        try:
            scan = find_disconnectable_bonds(product)
        except Exception:
            s.success += 1
            continue

        if not scan.get("ok"):
            s.success += 1
            continue

        precursors = None
        for bi in scan.get("bonds", []):
            tpls = bi.get("templates", [])
            if not tpls:
                continue
            br = execute_disconnection(product, tuple(bi["atoms"]), tpls[0]["name"])
            if br.success and br.precursors:
                precursors = br.precursors
                break

        if precursors is None:
            s.success += 1
            continue

        t0 = time.time()
        try:
            vf = validate_forward(precursors, product)
        except Exception as e:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": f"validate_forward exception: {e}", "row": rec.row_index})
            continue
        s.timings.append(time.time() - t0)

        if not vf.get("ok"):
            s.success += 1
            continue

        ok = True
        if "checks" not in vf or "assessment" not in vf:
            ok = False
            s.errors.append({"smiles": product, "reason": "missing checks/assessment", "row": rec.row_index})

        score = vf.get("assessment", {}).get("feasibility_score")
        if score is not None and not (0.0 <= score <= 1.0):
            ok = False
            s.errors.append({"smiles": product, "reason": f"score out of range: {score}", "row": rec.row_index})

        if ok:
            s.success += 1
        else:
            s.fail += 1


# ── Step 4: 端到端管线 ──────────────────────────────────────────
def run_pipeline(records, stats: Dict[str, Stats]):
    """完整管线: analyze → scan → break → validate"""
    s = stats["pipeline"]
    for rec in records:
        product = rec.product_smiles
        t0 = time.time()
        ok = True

        # analyze
        try:
            a = analyze_molecule(product)
            if "atoms" not in a:
                ok = False
                s.errors.append({"smiles": product, "reason": "analyze: missing atoms", "row": rec.row_index})
                s.fail += 1
                s.timings.append(time.time() - t0)
                continue
        except Exception as e:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": f"analyze exception: {e}", "row": rec.row_index})
            s.timings.append(time.time() - t0)
            continue

        # scan
        try:
            scan = find_disconnectable_bonds(product)
        except Exception as e:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": f"scan exception: {e}", "row": rec.row_index})
            s.timings.append(time.time() - t0)
            continue

        if not scan.get("ok") or not scan.get("bonds"):
            s.success += 1
            s.timings.append(time.time() - t0)
            continue

        # break (first bond with template)
        br = None
        for bi in scan["bonds"]:
            tpls = bi.get("templates", [])
            if tpls:
                try:
                    br = execute_disconnection(product, tuple(bi["atoms"]), tpls[0]["name"])
                except Exception:
                    pass
                break

        if br is None or not br.success or not br.precursors:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": "break failed or no precursors", "row": rec.row_index})
            s.timings.append(time.time() - t0)
            continue

        # validate
        try:
            vf = validate_forward(br.precursors, product)
            if vf.get("ok") and ("checks" not in vf or "assessment" not in vf):
                ok = False
                s.errors.append({"smiles": product, "reason": "validate: missing keys", "row": rec.row_index})
        except Exception as e:
            ok = False
            s.errors.append({"smiles": product, "reason": f"validate exception: {e}", "row": rec.row_index})

        s.timings.append(time.time() - t0)
        if ok:
            s.success += 1
        else:
            s.fail += 1


# ── Step 5: 逆合成还原评估 ──────────────────────────────────────
def run_retro(records, stats: Dict[str, Stats]):
    """穷举所有可断键位，比较生成前体与真实前体的 Tanimoto 相似度"""
    s = stats["retro"]
    tanimoto_scores: List[float] = []
    exact_matches = 0

    for rec in records:
        product = rec.product_smiles
        t0 = time.time()

        try:
            scan = find_disconnectable_bonds(product)
        except Exception:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": "scan exception", "row": rec.row_index})
            s.timings.append(time.time() - t0)
            continue

        if not scan.get("ok") or not scan.get("bonds"):
            s.fail += 1
            s.errors.append({"smiles": product, "reason": "no disconnectable bonds", **mol_features(product), "row": rec.row_index})
            s.timings.append(time.time() - t0)
            continue

        # 收集所有成功的前体集
        all_prec_sets = []
        for bi in scan["bonds"]:
            tpls = bi.get("templates", [])
            if not tpls:
                continue
            try:
                br = execute_disconnection(product, tuple(bi["atoms"]), tpls[0]["name"])
            except Exception:
                continue
            if br.success and br.precursors:
                all_prec_sets.append(br.precursors)

        if not all_prec_sets:
            s.fail += 1
            s.errors.append({"smiles": product, "reason": "all breaks failed", "row": rec.row_index})
            s.timings.append(time.time() - t0)
            continue

        # 真实前体
        true_frags = rec.precursors_smiles.split(".")
        true_mols = [m for f in true_frags if (m := parse_mol(f))]
        true_cans = {canonical(f) for f in true_frags if canonical(f)}

        best_sim = 0.0
        found_exact = False

        for prec_set in all_prec_sets:
            for gen_smi in prec_set:
                gen_mol = parse_mol(gen_smi)
                if gen_mol is None:
                    continue
                gc = canonical(gen_smi)
                if gc and gc in true_cans:
                    found_exact = True
                for tm in true_mols:
                    try:
                        sim = tanimoto(gen_mol, tm)
                        best_sim = max(best_sim, sim)
                    except Exception:
                        pass

        tanimoto_scores.append(best_sim)
        if found_exact:
            exact_matches += 1

        s.timings.append(time.time() - t0)
        s.success += 1

    # 附加统计到 stats
    if tanimoto_scores:
        tanimoto_scores.sort()
        n = len(tanimoto_scores)
        s.extra = {  # type: ignore[attr-defined]
            "tanimoto_mean": round(sum(tanimoto_scores) / n, 4),
            "tanimoto_median": round(tanimoto_scores[n // 2], 4),
            "tanimoto_p25": round(tanimoto_scores[n // 4], 4),
            "tanimoto_p75": round(tanimoto_scores[3 * n // 4], 4),
            "exact_match_count": exact_matches,
            "exact_match_rate": round(exact_matches / n, 4),
        }


# ── 报告输出 ─────────────────────────────────────────────────────
def print_report(file_label: str, stats: Dict[str, Stats], elapsed_total: float):
    print(f"\n{'='*60}")
    print(f"  验证报告: {file_label}")
    print(f"  总耗时: {elapsed_total:.1f}s")
    print(f"{'='*60}")

    for step_name in STEPS_ALL:
        if step_name not in stats:
            continue
        s = stats[step_name]
        if s.total == 0:
            continue
        status = "✓" if s.fail == 0 else "✗"
        print(f"\n  {status} {step_name:12s}  通过 {s.success}/{s.total} ({s.rate:.1%})  avg {s.avg_ms:.1f}ms")
        if s.errors:
            # 按 reason 聚合
            reason_counts: Dict[str, int] = defaultdict(int)
            for e in s.errors:
                reason_counts[e["reason"]] += 1
            top = sorted(reason_counts.items(), key=lambda x: -x[1])[:5]
            for reason, count in top:
                print(f"      [{count}x] {reason}")

        # retro 额外统计
        extra = getattr(s, "extra", None)
        if extra:
            print(f"      Tanimoto: mean={extra['tanimoto_mean']}, median={extra['tanimoto_median']}, "
                  f"P25={extra['tanimoto_p25']}, P75={extra['tanimoto_p75']}")
            print(f"      精确匹配: {extra['exact_match_count']}/{s.success} ({extra['exact_match_rate']:.1%})")

    # 失败详情（最多10条）
    all_errors = []
    for step_name, s in stats.items():
        for e in s.errors:
            all_errors.append({**e, "step": step_name})
    if all_errors:
        print(f"\n  失败案例 (前10条):")
        for e in all_errors[:10]:
            row = e.get("row", "?")
            print(f"    row-{row} [{e['step']}] {e['reason']}")
            if "smiles" in e:
                smi = e["smiles"]
                print(f"      SMILES: {smi[:80]}{'...' if len(smi) > 80 else ''}")

    print(f"\n{'='*60}\n")


# ── 主函数 ───────────────────────────────────────────────────────
STEP_RUNNERS = {
    "parse": run_parse,
    "bond": run_bond,
    "forward": run_forward,
    "pipeline": run_pipeline,
    "retro": run_retro,
}


def main():
    parser = argparse.ArgumentParser(description="化学层数据驱动快速验证")
    parser.add_argument("--file", default="BTF1",
                        help="数据文件简称 (BTF1-10, ATF, NEW, ALL)，默认 BTF1")
    parser.add_argument("--size", type=int, default=100,
                        help="每个文件采样数量，默认 100")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子，默认 42")
    parser.add_argument("--step", default=None,
                        help="要跑的步骤，逗号分隔 (parse,bond,forward,pipeline,retro)，默认全部")
    args = parser.parse_args()

    # 解析步骤
    if args.step:
        steps = [s.strip() for s in args.step.split(",")]
        for s in steps:
            if s not in STEP_RUNNERS:
                print(f"未知步骤: {s}，可选: {', '.join(STEP_RUNNERS.keys())}")
                sys.exit(1)
    else:
        steps = STEPS_ALL

    # 解析文件
    file_key = args.file.upper()
    if file_key == "ALL":
        files = list(FILE_MAP.items())
    elif file_key in FILE_MAP:
        files = [(file_key, FILE_MAP[file_key])]
    else:
        print(f"未知文件: {args.file}，可选: {', '.join(FILE_MAP.keys())}, ALL")
        sys.exit(1)

    print(f"化学层数据驱动验证")
    print(f"步骤: {', '.join(steps)}")
    print(f"文件: {', '.join(k for k, _ in files)}")
    print(f"采样: {args.size} 条/文件, seed={args.seed}")

    for file_label, file_path in files:
        if not file_path.exists():
            print(f"\n跳过 {file_label}: 文件不存在 ({file_path})")
            continue

        records = load_csv(str(file_path), sample_size=args.size, random_seed=args.seed)
        skipped = get_skipped_count()
        print(f"\n加载 {file_label}: {len(records)} 条记录 (跳过 {skipped} 空行)")

        stats: Dict[str, Stats] = {s: Stats() for s in steps}
        t_total = time.time()

        for step_name in steps:
            print(f"  运行 {step_name}...", end="", flush=True)
            t0 = time.time()
            STEP_RUNNERS[step_name](records, stats)
            print(f" {time.time() - t0:.1f}s")

        print_report(file_label, stats, time.time() - t_total)


if __name__ == "__main__":
    main()
