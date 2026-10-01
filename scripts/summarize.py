"""Aggregate run JSONs: CV mean/SD, fold vs seed variance, pooled confusion matrices, per-class P/R, LOSO per cohort,
Nadeau-Bengio comparisons. Output: results/summary/<config>.json"""
import json, glob, pathlib
import numpy as np, pandas as pd
from sklearn.metrics import precision_recall_fscore_support
from common import RES, write_json
from evalstats import nadeau_bengio, chance_levels

RUNS = RES / "runs"

def load(config, task, scheme="cv", tag=""):
    fs = sorted(glob.glob(str(RUNS / config / task / f"{scheme}_*.json")))
    out = []
    for f in fs:
        r = json.load(open(f))
        if r.get("tag", "") == tag: out.append(r)
    return out

def fold_table(config, task, metric="macro_f1", tag=""):
    rs = load(config, task, "cv", tag)
    return pd.DataFrame([{"fold": r["fold"], "seed": r["seed"], "v": r["metrics"].get(metric), "n_train": r["n_train"], "n_test": r["n_test"]} for r in rs])

def cv_summary(config, task, tag=""):
    rs = load(config, task, "cv", tag)
    if not rs: return None
    df = pd.DataFrame([{"fold": r["fold"], "seed": r["seed"], "f1": r["metrics"]["macro_f1"], "mcc": r["metrics"]["mcc"],
                        "auc": r["metrics"].get("auc")} for r in rs])
    per_fold = df.groupby("fold")[["f1", "mcc", "auc"]].mean()
    seed_sd = df.groupby("fold").f1.std(ddof=1).mean() if df.seed.nunique() > 1 else 0.0
    labels = rs[0]["metrics"]["confusion_matrix"]["labels"]
    cm = np.sum([np.array(r["metrics"]["confusion_matrix"]["matrix"]) for r in rs], 0)
    y = np.concatenate([r["y"] for r in rs]); p = np.concatenate([r["pred"] for r in rs])
    P, R, F, N = precision_recall_fscore_support(y, p, labels=labels, zero_division=0)
    return {"n_folds": int(per_fold.shape[0]), "n_seeds": int(df.seed.nunique()),
            "macro_f1_mean": float(per_fold.f1.mean()), "macro_f1_fold_sd": float(per_fold.f1.std(ddof=1)),
            "macro_f1_seed_sd": float(seed_sd), "mcc_mean": float(per_fold.mcc.mean()), "mcc_fold_sd": float(per_fold.mcc.std(ddof=1)),
            "auc_mean": float(per_fold.auc.mean()), "auc_fold_sd": float(per_fold.auc.std(ddof=1)),
            "per_class_pooled": {l: {"precision": float(P[i]), "recall": float(R[i]), "f1": float(F[i])} for i, l in enumerate(labels)},
            "confusion_matrix_summed_over_folds_and_seeds": {"labels": labels, "matrix": cm.tolist()},
            "fold_sizes": {"n_train_mean": float(np.mean([r["n_train"] for r in rs])), "n_test_mean": float(np.mean([r["n_test"] for r in rs]))},
            "per_fold_macro_f1": per_fold.f1.round(4).to_dict()}

def loso_summary(config, task, tag=""):
    rs = load(config, task, "loso", tag)
    out = {}
    for c in sorted({r["fold"] for r in rs}):
        rr = [r for r in rs if r["fold"] == c]
        f1 = [r["metrics"]["macro_f1"] for r in rr]
        m = rr[0]["metrics"]
        ch = chance_levels(rr[0]["y"], n_sim=1000)
        out[c.replace("loso_", "")] = {"n_test": rr[0]["n_test"], "macro_f1": float(np.mean(f1)), "macro_f1_seed_sd": float(np.std(f1, ddof=1)) if len(f1) > 1 else 0.0,
                                       "mcc": m["mcc"], "auc": m.get("auc"), "per_class": m["per_class"], "confusion_matrix": m["confusion_matrix"],
                                       "chance_uniform_macroF1": round(ch["uniform_random_macroF1"], 3), "chance_majority_macroF1": round(ch["majority_class_macroF1"], 3)}
    return out

def compare(cfg_a, cfg_b, task, metric="macro_f1", tag_a="", tag_b=""):
    """Nadeau-Bengio on the 25 fold-level differences (seed-averaged per fold), a minus b."""
    A = fold_table(cfg_a, task, metric, tag_a).groupby("fold").agg(v=("v", "mean"), ntr=("n_train", "mean"), nte=("n_test", "mean"))
    B = fold_table(cfg_b, task, metric, tag_b).groupby("fold").v.mean()
    j = A.join(B.rename("w"), how="inner")
    if len(j) < 2: return None
    return nadeau_bengio((j.v - j.w).values, j.ntr.mean(), j.nte.mean())

def summarize(config, tasks, tag=""):
    out = {"config": config, "tag": tag, "tasks": {}}
    ch = json.load(open(RES / "dev_task_sizes_chance.json"))
    for t in tasks:
        out["tasks"][t] = {"cv": cv_summary(config, t, tag), "loso": loso_summary(config, t, tag),
                           "chance_dev": {k: ch[t][k] for k in ("uniform_random_macroF1", "prior_random_macroF1", "majority_class_macroF1")}}
    rs = glob.glob(str(RUNS / config / "*" / "*.json"))
    if rs: out["versions"] = json.load(open(rs[0]))["versions"]
    write_json(out, RES / "summary" / f"{config}{tag}.json")
    return out

def brief(out):
    rows = []
    for t, v in out["tasks"].items():
        cv = v["cv"]
        if not cv: continue
        lo = v["loso"]
        rows.append({"task": t, "cv_macroF1": f'{cv["macro_f1_mean"]:.3f} ± {cv["macro_f1_fold_sd"]:.3f}', "seed_sd": round(cv["macro_f1_seed_sd"], 3),
                     "cv_MCC": round(cv["mcc_mean"], 3), "cv_AUC": round(cv["auc_mean"], 3),
                     "chance_uniform": round(v["chance_dev"]["uniform_random_macroF1"], 3),
                     "LOSO_macroF1": {k: round(x["macro_f1"], 3) for k, x in lo.items()}})
    return pd.DataFrame(rows)
