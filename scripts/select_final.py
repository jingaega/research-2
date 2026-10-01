"""Final selection per PLAN §9 + bootstrap bias-corrected CV (BBC-CV, Tsamardinos 2018) over the candidate pool.
Writes results/variants.json (full variant table with running count) and results/final_selection.json."""
import json, glob, numpy as np, pandas as pd
from summarize import cv_summary
from pipeline import TASK_ORDER, RUNS
from sklearn.metrics import f1_score
from common import RES, write_json
VARIANTS = [  # (id, config, counts_against_budget, family, eligible_final, verdict)
 ("C0", "C0_baseline_LR", False, "simple", True, "control: no-structure baseline"),
 ("C1", "C1_study_only", False, "control", False, "control: study-only"),
 ("V5", "V5_combat_LR", True, "batch-correction", False, "INVALID: composition leakage (C8 within-batch permutation)"),
 ("V1", "V1_string", True, "graph-curated", True, "FAILED: no gain vs no-graph or random graph"),
 ("V2", "V2_coexpr", True, "graph-expression", True, "FAILED: no gain vs no-graph or random graph"),
 ("V3", "V3_reactome", True, "pathway-curated", True, "FAILED: worse than no-graph; random sets better on 3/4 tasks"),
 ("V6a", "V6a_filter1000_LR", True, "simple", True, "no gain beyond noise (+0.014..+0.023, NS)"),
 ("V6b", "V6b_filter1000_STRING", False, "graph-curated", True, "equal-treatment twin of V6 (counted once with V6)"),
 ("V7", "V7_RF", True, "simple", True, "FAILED: worse on BD/MDD, NS gains SCZ/B"),
 ("V8", "V8_ensemble_C0_V1", True, "ensemble", True, "FAILED: = C0"),
 ("V9", "V9_celltype", True, "cell-type composition", True, "FAILED as model; recovers 88% of C0 SCZ excess over chance"),
 ("V4", None, False, "GNN", False, "NOT RUN: pre-registered trigger not met"),
]
rows, count = [], 0
for vid, cfg, counts, fam, elig, verdict in VARIANTS:
    if counts: count += 1
    r = {"id": vid, "config": cfg, "family": fam, "counts_against_budget": counts, "running_count": count if counts else None,
         "eligible_for_final": elig, "verdict": verdict}
    if cfg:
        for t in TASK_ORDER:
            s = cv_summary(cfg, t); r[t] = round(s["macro_f1_mean"], 4); r[t + "_sd"] = round(s["macro_f1_fold_sd"], 4)
        r["mean_4_tasks"] = round(np.mean([r[t] for t in TASK_ORDER]), 4)
    rows.append(r)
V = pd.DataFrame(rows)
print(V[["id", "family", "running_count", "eligible_for_final", "mean_4_tasks"] + TASK_ORDER].to_string())
el = V[V.eligible_for_final & V.config.notna()]
simple = el[el.id.isin(["C0", "V6a", "V7"])]
nonbase = el[~el.id.isin(["C0", "V6a", "V7"])]
best_simple = simple.sort_values("mean_4_tasks", ascending=False).iloc[0]
best_model = nonbase.sort_values("mean_4_tasks", ascending=False).iloc[0]
print("best simple baseline:", best_simple.id, best_simple.mean_4_tasks, "| best model:", best_model.id, best_model.mean_4_tasks)

# ---- BBC-CV over the eligible pool (selection by mean macro-F1 over tasks), per repeat, 1000 bootstraps
def oof(cfg, t, rep):
    P = {}
    for f in glob.glob(str(RUNS / cfg / t / f"cv_r{rep}f*_s0.json")):
        r = json.load(open(f))
        for d, y, p in zip(r["donors"], r["y"], r["pred"]): P[d] = (y, p)
    return P
pool = el.config.tolist(); rng = np.random.default_rng(0); bbc = {c: [] for c in ["selected_mean"] + TASK_ORDER}
for rep in range(5):
    O = {t: {c: oof(c, t, rep) for c in pool} for t in TASK_ORDER}
    ids = {t: sorted(O[t][pool[0]]) for t in TASK_ORDER}
    Y = {t: np.array([O[t][pool[0]][d][0] for d in ids[t]]) for t in TASK_ORDER}
    Pr = {t: {c: np.array([O[t][c][d][1] for d in ids[t]]) for c in pool} for t in TASK_ORDER}
    for b in range(200):
        bi = {t: rng.integers(0, len(ids[t]), len(ids[t])) for t in TASK_ORDER}
        oob = {t: np.setdiff1d(np.arange(len(ids[t])), bi[t]) for t in TASK_ORDER}
        sc = {c: np.mean([f1_score(Y[t][bi[t]], Pr[t][c][bi[t]], average="macro") for t in TASK_ORDER]) for c in pool}
        best = max(sc, key=sc.get)
        per = [f1_score(Y[t][oob[t]], Pr[t][best][oob[t]], average="macro") for t in TASK_ORDER]
        bbc["selected_mean"].append(np.mean(per))
        for t, v in zip(TASK_ORDER, per): bbc[t].append(v)
bbc_res = {k: float(np.mean(v)) for k, v in bbc.items()}
print("BBC-CV estimate of the selected configuration (mean over tasks):", round(bbc_res["selected_mean"], 4), {t: round(bbc_res[t], 3) for t in TASK_ORDER})
write_json({"variants": rows, "running_count_final": count, "budget": 12}, RES / "variants.json")
write_json({"rule": "PLAN §9: best model = highest mean dev-CV macro-F1 over 4 tasks among eligible non-baseline configs; best simple baseline = best of {C0, V6a, V7}",
            "best_model": {"id": best_model.id, "config": best_model.config, "mean_cv": best_model.mean_4_tasks},
            "best_simple_baseline": {"id": best_simple.id, "config": best_simple.config, "mean_cv": best_simple.mean_4_tasks},
            "excluded": {"V5": "composition leakage (C8)"},
            "bbc_cv_pool": pool, "bbc_cv_estimate_of_selected": bbc_res,
            "note": "BBC-CV: bootstrap donors' pooled out-of-fold predictions per repeat, select the best config on the bootstrap sample, evaluate on out-of-bag donors (200 boots x 5 repeats)."},
           RES / "final_selection.json")
