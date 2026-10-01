"""Controls C4(i) sample-level vs donor-grouped split, and C6 sex-mismatch sensitivity (C0 model).

C4(i): multi-sample-per-donor SCZ-vs-CTL dataset built from DEV donors only (never test donors): every usable array
of the SCZ-task cohorts (Stanley parietal + cerebellum + GSE12649 BA46, Pittsburgh PFC/STR/BA9/BA25, Pritzker per array,
CharingCross, Victoria). Batch unit = series|chip|region. Same C0 model; (a) StratifiedKFold on arrays (ignores donors)
vs (b) StratifiedGroupKFold by donor split group. 5 x 5 each.
C6: base = all eligible dev donors (incl. affected); 'minus affected' = primary C0; 'minus random' = drop donors
matched to the affected ones on diagnosis x cohort, 20 draws. 5x5 CV each."""
import glob, json, pathlib
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from pipeline import run_config, run_fold, batch_std, fit_predict, load_dev, task_data, TASK_ORDER, TASKS
from summarize import summarize, brief, compare, cv_summary
from evalstats import metrics
from common import DATA, RES, write_json
pd.set_option("display.width", 250)

# ------------------------------------------------------------- C4(i)
X, meta = load_dev()
smap = pd.read_csv(DATA / "sample_donor_map.tsv", sep="\t", dtype=str, keep_default_na=False)
genes = list(X.columns)
dev_ok = meta[~meta.affected & meta.dx.isin(["SCZ", "CTL"]) & meta.cohort.isin(TASKS["A_SCZ_vs_CTL"]["cohorts"])]
s = smap[smap.donor_id.isin(dev_ok.index) & ~smap.region.isin(["HPC"])].copy()
arrays = {}
for f in glob.glob(str(DATA / "expr" / "*.pkl")):
    g = pd.read_pickle(f); name = pathlib.Path(f).stem
    for c in g.columns:
        if c in set(s.gsm): arrays[c] = (g[c].reindex(genes).values, name)
s = s[s.gsm.isin(arrays)]
s["batch"] = [arrays[g][1] + "|" + r for g, r in zip(s.gsm, s.region)]
s = s[s.groupby("batch").gsm.transform("size") >= 10]          # need enough arrays per batch for standardisation
XA = np.stack([arrays[g][0] for g in s.gsm]); assert not np.isnan(XA).any()
s["dx_d"] = dev_ok.loc[s.donor_id, "dx"].values; s["grp"] = dev_ok.loc[s.donor_id, "split_group"].values
print("C4(i) arrays:", len(s), "donors:", s.donor_id.nunique(), "batches:", s.batch.nunique())
y = s.dx_d.values; strat = y + "|" + s.batch.values
out = {}
for scheme in ["sample_split_LEAK", "donor_grouped"]:
    f1s, aucs = [], []
    for r in range(5):
        spl = StratifiedKFold(5, shuffle=True, random_state=r).split(XA, strat) if scheme.startswith("sample") else \
              StratifiedGroupKFold(5, shuffle=True, random_state=r).split(XA, strat, s.grp)
        for tr, te in spl:
            Xtr, Xte = batch_std(XA[tr], s.batch.values[tr], XA[te], s.batch.values[te])
            cl, P, _ = fit_predict("lr", Xtr, y[tr], Xte, s.grp.values[tr], r)
            m = metrics(y[te], cl[P.argmax(1)], P, cl); f1s.append(m["macro_f1"]); aucs.append(m["auc"])
            if scheme.startswith("sample"):
                shared = len(set(s.donor_id.values[tr]) & set(s.donor_id.values[te]))
                out.setdefault("donors_shared_train_test_per_fold", []).append(shared)
    out[scheme] = {"macro_f1_mean": float(np.mean(f1s)), "macro_f1_fold_sd": float(np.std(f1s, ddof=1)), "auc_mean": float(np.mean(aucs)), "per_fold": f1s}
    print("C4(i)", scheme, {k: v for k, v in out[scheme].items() if k != "per_fold"})
from evalstats import nadeau_bengio
out["nadeau_bengio_sample_minus_grouped"] = nadeau_bengio(np.array(out["sample_split_LEAK"]["per_fold"]) - np.array(out["donor_grouped"]["per_fold"]), 0.8 * len(s), 0.2 * len(s))
out["note"] = "per-fold pairing is by index only (different splits); the NB p is indicative. n_arrays=%d, n_donors=%d" % (len(s), s.donor_id.nunique())
out["n_arrays"] = int(len(s)); out["n_donors"] = int(s.donor_id.nunique()); out["arrays_per_donor"] = s.groupby("donor_id").size().value_counts().sort_index().to_dict()
write_json(out, RES / "summary" / "C4a_sample_vs_donor_split.json")

# ------------------------------------------------------------- C6
spec = {"model": "lr", "adjust": "bstd"}
run_config("C6_all_eligible", spec, population="all")
print("C6 all eligible"); print(brief(summarize("C6_all_eligible", TASK_ORDER)).to_string())
rng = np.random.default_rng(2026)
c6 = {}
for t in TASK_ORDER:
    Xa, ma = task_data(t, population="all")
    aff = ma[ma.affected]
    base = cv_summary("C6_all_eligible", t)["macro_f1_mean"]; prim = cv_summary("C0_baseline_LR", t)["macro_f1_mean"]
    c6[t] = {"n_affected_in_task": int(len(aff)), "affected_composition": aff.groupby(["dx", "cohort"]).size().rename(lambda x: "|".join(x)).to_dict(),
             "all_eligible": base, "minus_affected(primary)": prim, "random_drops": []}
    if len(aff) == 0: continue
    for k in range(20):
        drop = []
        for (dx, coh), n in aff.groupby(["dx", "cohort"]).size().items():
            pool = ma[(ma.dx == dx) & (ma.cohort == coh) & ~ma.affected].index
            drop += list(rng.choice(pool, n, replace=False))
        cfg = f"C6_random_drop"
        run_config(cfg, spec, tasks=[t], schemes=("cv",), population="all", drop_ids=drop, tag=f"_d{k}")
        c6[t]["random_drops"].append(cv_summary(cfg, t, tag=f"_d{k}")["macro_f1_mean"])
    rd = np.array(c6[t]["random_drops"])
    c6[t].update({"random_drop_mean": float(rd.mean()), "random_drop_sd": float(rd.std(ddof=1)),
                  "delta_minus_affected": prim - base, "delta_random_mean": float(rd.mean() - base),
                  "frac_random_drops_with_delta_ge_affected": float(np.mean(rd - base >= prim - base))})
    print("C6", t, {k: v for k, v in c6[t].items() if k != "random_drops"})
write_json(c6, RES / "summary" / "C6_sex_mismatch_sensitivity.json")
