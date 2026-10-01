"""C9 sensitivity: C0 and V5 on each task restricted to batches that contain every class of the task."""
import pandas as pd
from pipeline import run_config, task_data, TASK_ORDER
from summarize import summarize, brief, compare
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
for t in TASK_ORDER:
    X, m = task_data(t)
    k = m.groupby("batch").dx.nunique(); bad = k[k < m.dx.nunique()].index
    drop = list(m.index[m.batch.isin(bad)])
    print(t, "dropping", len(drop), "donors from batches", list(bad))
    for cfg, spec in [("C9_C0_completebatches", {"model": "lr", "adjust": "bstd"}), ("C9_V5_completebatches", {"model": "lr", "adjust": "combat_infold"})]:
        run_config(cfg, spec, tasks=[t], drop_ids=drop)
for cfg in ["C9_C0_completebatches", "C9_V5_completebatches"]:
    print(cfg); print(brief(summarize(cfg, TASK_ORDER)).to_string())
for t in TASK_ORDER:
    r = compare("C9_V5_completebatches", "C9_C0_completebatches", t)
    print(t, "V5-C0 on complete batches: Δ %+.3f p_corr %.3f p_unc %.4f" % (r["mean_diff"], r["p_corrected"], r["p_uncorrected"]))
