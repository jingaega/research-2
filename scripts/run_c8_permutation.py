"""C8: within-BATCH label permutation null (preserves batch x class counts) for C0 and V5 (all 4 tasks, 5 permutations x 5x5 CV)."""
import numpy as np, pandas as pd, json
import pipeline as P
from summarize import cv_summary
from common import RES, write_json
orig = P.task_data
out = {}
for k in range(5):
    def permuted(task, population="primary", drop_ids=(), exclude_cohorts=(), _k=k):
        X, meta = orig(task, population, drop_ids, exclude_cohorts)
        meta = meta.copy(); rng = np.random.default_rng(1000 + _k)
        for c in meta.batch.unique():
            i = np.where(meta.batch.values == c)[0]
            meta.iloc[i, meta.columns.get_loc("dx")] = rng.permutation(meta.dx.values[i])
        return X, meta
    P.task_data = permuted
    for cfg, spec in [("C8_perm_C0", {"model": "lr", "adjust": "bstd"}), ("C8_perm_V5", {"model": "lr", "adjust": "combat_infold"})]:
        P.run_config(cfg, spec, schemes=("cv",), tag=f"_p{k}")
P.task_data = orig
for cfg in ["C8_perm_C0", "C8_perm_V5"]:
    out[cfg] = {}
    for t in P.TASK_ORDER:
        v = [cv_summary(cfg, t, tag=f"_p{k}")["macro_f1_mean"] for k in range(5)]
        a = [cv_summary(cfg, t, tag=f"_p{k}")["auc_mean"] for k in range(5)]
        out[cfg][t] = {"macro_f1_per_perm": v, "macro_f1_mean": float(np.mean(v)), "auc_mean": float(np.mean(a))}
        print(cfg, t, "macroF1 %.3f (perms %s) AUC %.3f" % (np.mean(v), np.round(v, 3), np.mean(a)))
write_json(out, RES / "summary" / "C8_within_batch_permutation.json")
