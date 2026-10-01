"""V9 cell-type composition classifier (+ random gene-set control, 3 seeds)."""
import json, pandas as pd
from pipeline import run_config, TASK_ORDER
from summarize import summarize, brief, compare, cv_summary
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
run_config("V9_celltype", {"model": "lr", "adjust": "bstd", "transform": "celltype"})
run_config("V9_celltype_RANDOM", {"model": "lr", "adjust": "bstd", "transform": "celltype_random"}, seeds=(0, 1, 2))
for c in ["V9_celltype", "V9_celltype_RANDOM"]:
    print(c); print(brief(summarize(c, TASK_ORDER)).to_string())
ch = json.load(open("../results/dev_task_sizes_chance.json"))
for t in TASK_ORDER:
    a = compare("V9_celltype", "C0_baseline_LR", t); b = compare("V9_celltype", "V9_celltype_RANDOM", t)
    v9 = cv_summary("V9_celltype", t)["macro_f1_mean"]; c0 = cv_summary("C0_baseline_LR", t)["macro_f1_mean"]; u = ch[t]["uniform_random_macroF1"]
    print(f"{t}: V9-C0 Δ {a['mean_diff']:+.3f} p_corr {a['p_corrected']:.3f} | V9-random Δ {b['mean_diff']:+.3f} p_corr {b['p_corrected']:.3f} | "
          f"fraction of C0 excess over chance: {(v9-u)/(c0-u):.2f}")
