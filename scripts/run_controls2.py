"""Controls: C4(ii) ComBat leakage demo (+ V5 in-fold arm), C7 single-study dependence, C2 study-identity probe."""
import json, numpy as np, pandas as pd
from pipeline import run_config, task_data, cv_folds, batch_std, combat, fit_predict, load_dev, TASK_ORDER
from summarize import summarize, brief, compare
from evalstats import metrics, chance_levels
from common import RES, write_json
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)

# ---- V5 / C4(ii)
run_config("V5_combat_LR", {"model": "lr", "adjust": "combat_infold"})
run_config("C4b_combat_alldata_LEAK", {"model": "lr", "adjust": "combat_alldata"}, schemes=("cv",))
for c in ["V5_combat_LR", "C4b_combat_alldata_LEAK"]:
    print(c); print(brief(summarize(c, TASK_ORDER)).to_string())
for t in TASK_ORDER:
    for a, b in [("V5_combat_LR", "C0_baseline_LR"), ("C4b_combat_alldata_LEAK", "V5_combat_LR")]:
        r = compare(a, b, t); print(t, a, "vs", b, "diff %.3f p_corr %.3f p_unc %.4f" % (r["mean_diff"], r["p_corrected"], r["p_uncorrected"]))

# ---- C7: BD single-study dependence (Stanley supplies >= 50% of dev BD donors)
X, meta = load_dev(); bd = meta[(meta.dx == "BD") & ~meta.affected]
share = bd.cohort.value_counts(normalize=True); print("dev BD share by cohort:", share.round(3).to_dict())
dom = share.index[share >= 0.5].tolist()
for c in dom:
    run_config(f"C7_C0_without_{c}", {"model": "lr", "adjust": "bstd"}, tasks=["A_BD_vs_CTL", "B_multiclass"], exclude_cohorts=(c,))
    print("C7 without", c); print(brief(summarize(f"C7_C0_without_{c}", ["A_BD_vs_CTL", "B_multiclass"])).to_string())

# ---- C2: study-identity probe (predict cohort / batch from features), union of all primary dev task donors
m = ~meta.affected & meta.dx.isin(["SCZ", "BD", "MDD", "CTL"])
Xp, mp = X[m], meta[m]
res = {}
from sklearn.model_selection import StratifiedGroupKFold
for target in ["cohort", "batch"]:
    y = mp[target].values
    for adj in ["none", "bstd", "combat"]:
        f1s = []
        for r in range(5):
            for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=r).split(Xp, y, mp.split_group):
                Xtr, Xte = Xp.values[tr], Xp.values[te]
                if adj == "none":
                    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6; Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
                elif adj == "bstd": Xtr, Xte = batch_std(Xtr, mp.batch.values[tr], Xte, mp.batch.values[te])
                else: Xtr, Xte = combat(Xtr, mp.batch.values[tr], mp.dx.values[tr], Xte, mp.batch.values[te])
                cl, P, _ = fit_predict("lr", Xtr, y[tr], Xte, mp.split_group.values[tr], r)
                f1s.append(metrics(y[te], cl[P.argmax(1)], P, cl)["macro_f1"])
        ch = chance_levels(y, n_sim=2000)
        res[f"{target}|{adj}"] = {"macro_f1_mean": float(np.mean(f1s)), "macro_f1_fold_sd": float(np.std(f1s, ddof=1)),
                                  "chance_uniform": ch["uniform_random_macroF1"], "chance_majority": ch["majority_class_macroF1"], "n_classes": int(len(set(y)))}
        print("probe", target, adj, res[f"{target}|{adj}"])
write_json({"n_donors": int(len(Xp)), "results": res, "note": "adjustment fitted in-fold; test donors adjusted with their own (known) batch label, as in the pipeline"},
           RES / "summary" / "C2_study_identity_probe.json")
