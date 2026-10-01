"""V6 (ANOVA filter on C0 and V1), V7 (RF), V8 (ensemble C0+V1 from stored OOF probabilities)."""
import json, glob, numpy as np, pandas as pd
from pipeline import run_config, TASK_ORDER, RUNS
from summarize import summarize, brief, compare
from evalstats import metrics
from common import write_json
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
run_config("V6a_filter1000_LR", {"model": "lr", "adjust": "bstd", "filter_k": 1000})
run_config("V6b_filter1000_STRING", {"model": "lr", "adjust": "bstd", "transform": "string", "filter_k": 1000})
run_config("V7_RF", {"model": "rf", "adjust": "bstd"}, seeds=(0, 1, 2))
# V8: average stored probabilities fold by fold
for t in TASK_ORDER:
    for f in sorted(glob.glob(str(RUNS / "C0_baseline_LR" / t / "*.json"))):
        a = json.load(open(f)); b = json.load(open(f.replace("C0_baseline_LR", "V1_string")))
        assert a["donors"] == b["donors"] and a["classes"] == b["classes"]
        P = (np.array(a["proba"]) + np.array(b["proba"])) / 2; cl = np.array(a["classes"])
        r = dict(a); r.update({"config": "V8_ensemble_C0_V1", "spec": {"ensemble": ["C0_baseline_LR", "V1_string"]},
                               "proba": np.round(P, 5).tolist(), "pred": list(cl[P.argmax(1)]), "metrics": metrics(a["y"], cl[P.argmax(1)], P, cl)})
        write_json(r, f.replace("C0_baseline_LR", "V8_ensemble_C0_V1"))
for c in ["V6a_filter1000_LR", "V6b_filter1000_STRING", "V7_RF", "V8_ensemble_C0_V1"]:
    print(c); print(brief(summarize(c, TASK_ORDER)).to_string(), flush=True)
for t in TASK_ORDER:
    for a, b in [("V6a_filter1000_LR", "C0_baseline_LR"), ("V6b_filter1000_STRING", "V1_string"), ("V7_RF", "C0_baseline_LR"), ("V8_ensemble_C0_V1", "C0_baseline_LR")]:
        r = compare(a, b, t); print(f"{t} {a} vs {b}: Δ {r['mean_diff']:+.3f} p_corr {r['p_corrected']:.3f} (unc {r['p_uncorrected']:.4f})", flush=True)
