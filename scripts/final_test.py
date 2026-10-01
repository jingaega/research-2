"""SINGLE final evaluation on the locked test set (rule 11). Reads data/features_test.pkl for the first time."""
import json, numpy as np, pandas as pd
import pipeline as P
from evalstats import metrics, bootstrap_ci, chance_levels
from common import DATA, RES, write_json
from sklearn.metrics import f1_score, roc_auc_score, matthews_corrcoef
FINAL = {"V6b_best_model": {"model": "lr", "adjust": "bstd", "transform": "string", "filter_k": 1000},
         "V6a_best_simple_baseline": {"model": "lr", "adjust": "bstd", "filter_k": 1000},
         "C1_study_only_reference": {"model": "study_only"}}
Xd, md = P.load_dev()
Xt = pd.read_pickle(DATA / "features_test.pkl")
mt = pd.read_csv(DATA / "split_locked.tsv", sep="\t", index_col=0, keep_default_na=False).loc[Xt.index]
mt["affected"] = mt.affected.astype(str) == "True"
assert set(Xt.index).isdisjoint(Xd.index) and (Xt.columns == Xd.columns).all()
out = {"versions": P.vers(), "protocol": "LOG.md 2026-10-01 frozen test protocol", "tasks": {}}
for t in P.TASK_ORDER:
    spec = P.TASKS[t]
    dtr = md[md.cohort.isin(spec["cohorts"]) & md.dx.isin(spec["classes"]) & ~md.affected]
    dte = mt[mt.cohort.isin(spec["cohorts"]) & mt.dx.isin(spec["classes"]) & ~mt.affected]
    assert set(dte.batch) <= set(dtr.batch), "unseen batch in test"
    X = pd.concat([Xd.loc[dtr.index], Xt.loc[dte.index]]); M = pd.concat([dtr, dte])
    tr = np.arange(len(dtr)); te = np.arange(len(dtr), len(M))
    k = M.iloc[tr].groupby("batch").dx.nunique(); complete = set(k[k == M.dx.nunique()].index)
    sub = M.iloc[te].batch.isin(complete).values
    res_t = {"n_train": int(len(tr)), "n_test": int(len(te)), "test_class_n": dte.dx.value_counts().to_dict(),
             "test_cohort_dx": {f"{c}|{d}": int(n) for (c, d), n in dte.groupby(["cohort", "dx"]).size().items()},
             "chance_test": chance_levels(dte.dx.values, dtr.dx.values, n_sim=5000), "models": {}}
    preds = {}
    for name, sp in FINAL.items():
        r = P.run_fold("FINAL_TEST_" + name, sp, t, "test", "locked", tr, te, X, M, 0)
        y, p, Pr, cl = np.array(r["y"]), np.array(r["pred"]), np.array(r["proba"]), r["classes"]
        preds[name] = (p, Pr, cl)
        m = r["metrics"]; ci = bootstrap_ci(y, p, Pr, n_boot=2000, seed=0, classes=np.array(cl))
        msub = metrics(y[sub], p[sub], Pr[sub], np.array(cl)) if sub.sum() > 5 and len(set(y[sub])) == len(cl) else None
        res_t["models"][name] = {"metrics": m, "ci95": ci, "C": r["info"].get("C"),
                                 "complete_batch_subset": {"n": int(sub.sum()), "metrics": msub}}
    # paired bootstrap graph effect V6b - V6a
    y = dte.dx.values; rng = np.random.default_rng(1); d = []
    pa, pb = preds["V6b_best_model"][0], preds["V6a_best_simple_baseline"][0]
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) < len(set(y)): continue
        d.append(f1_score(y[i], pa[i], average="macro") - f1_score(y[i], pb[i], average="macro"))
    res_t["graph_effect_V6b_minus_V6a"] = {"macro_f1_diff": float(f1_score(y, pa, average="macro") - f1_score(y, pb, average="macro")),
                                           "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
                                           "frac_boot_le_0": float(np.mean(np.array(d) <= 0))}
    out["tasks"][t] = res_t
    print(t, "n_test", len(te), {n: round(v["metrics"]["macro_f1"], 3) for n, v in res_t["models"].items()},
          "chance(uniform)", round(res_t["chance_test"]["uniform_random_macroF1"], 3))
write_json(out, RES / "final_test_results.json")
