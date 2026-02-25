"""临时脚本: 展示 0.5-0.7 区间的案例"""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from Rachel.tools._near_miss import evaluate_single, FILE_MAP
from Rachel.tests.data_driven.loader import load_csv

records = load_csv(str(FILE_MAP["BTF1"]), sample_size=200, random_seed=42)
print(f"加载 {len(records)} 条, 评估中...")
cases = [evaluate_single(rec) for rec in records]

mid = [c for c in cases if 0.5 <= c.combined_best_tanimoto < 0.7]
mid.sort(key=lambda c: c.combined_best_tanimoto)

print(f"\n0.5-0.7 区间: {len(mid)} 条\n")

for i, c in enumerate(mid):
    fgi_tag = ""
    if c.best_fgi_tanimoto > c.best_tanimoto:
        fgi_tag = f"  [FGI提升: {c.best_tanimoto:.2f}->{c.best_fgi_tanimoto:.2f}]"
    print(f"{i+1:2d}. row={c.row_index:5d}  t={c.combined_best_tanimoto:.4f}  "
          f"rings={c.n_rings} atoms={c.n_heavy} [{c.heteroatoms}]{fgi_tag}")
    print(f"    产物: {c.product[:85]}")
    print(f"    真实: {c.true_precursors[:85]}")
    if c.best_gen_precursors:
        gen = ".".join(c.best_gen_precursors)[:85]
        print(f"    断键: {gen}  (t={c.best_tanimoto})")
        print(f"      tpl: {c.best_template[:65]}")
    if c.best_fgi_precursors:
        fgi = ".".join(c.best_fgi_precursors)[:85]
        print(f"    FGI:  {fgi}  (t={c.best_fgi_tanimoto})")
        print(f"      tpl: {c.best_fgi_template}")
    print()
