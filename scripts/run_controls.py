"""Required controls C0 (no-structure baseline) and C1 (study-only)."""
from pipeline import run_config
from summarize import summarize, brief
import pandas as pd
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
for cfg, spec in [("C0_baseline_LR", {"model": "lr", "adjust": "bstd"}), ("C1_study_only", {"model": "study_only"})]:
    n = run_config(cfg, spec)
    print(cfg, "fold-runs:", n); print(brief(summarize(cfg, ["A_SCZ_vs_CTL", "A_BD_vs_CTL", "A_MDD_vs_CTL", "B_multiclass"])).to_string())
