"""V1 (STRING), V2 (co-expression kNN), V3 (Reactome) with density-matched random controls (3 seeds).
Base adjustment = per-batch standardisation (as pre-registered; no-graph ablation = C0)."""
import sys
from pipeline import run_config, TASK_ORDER
from summarize import summarize, brief, compare
import pandas as pd
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
adj = sys.argv[1] if len(sys.argv) > 1 else "bstd"
suffix = "" if adj == "bstd" else "_combat"
base = "C0_baseline_LR" if adj == "bstd" else "V5_combat_LR"
plan = [(f"V1_string{suffix}", "string"), (f"V1_string_RANDOM{suffix}", "string_random"),
        (f"V3_reactome{suffix}", "reactome"), (f"V3_reactome_RANDOM{suffix}", "reactome_random"),
        (f"V2_coexpr{suffix}", "coexpr"), (f"V2_coexpr_RANDOM{suffix}", "coexpr_random")]
for cfg, tf in plan:
    seeds = (0, 1, 2) if tf.endswith("random") else (0,)
    run_config(cfg, {"model": "lr", "adjust": adj, "transform": tf}, seeds=seeds)
    print(cfg); print(brief(summarize(cfg, TASK_ORDER)).to_string(), flush=True)
for v in ["V1_string", "V3_reactome", "V2_coexpr"]:
    for t in TASK_ORDER:
        a = compare(v + suffix, base, t); b = compare(v + suffix, v + "_RANDOM" + suffix, t)
        print(f"{v}{suffix} {t}: vs no-graph Δ {a['mean_diff']:+.3f} p_corr {a['p_corrected']:.3f} (unc {a['p_uncorrected']:.3f}) | "
              f"vs random Δ {b['mean_diff']:+.3f} p_corr {b['p_corrected']:.3f} (unc {b['p_uncorrected']:.3f})", flush=True)
