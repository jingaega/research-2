"""Lock the final test set BEFORE any modeling: ~25% of eligible donors, grouped by split_group (liberal
duplicate groups), stratified by diagnosis x cohort. Uses all eligible donors, including 'affected'
ones, so the same split serves the primary analysis (affected removed) and the sensitivity analysis.
Output: data/split_locked.tsv, results/test_split.json"""
import json
import numpy as np, pandas as pd, sklearn
from sklearn.model_selection import StratifiedGroupKFold
from common import DATA, RES, write_json, versions

SEED = 20260930
D = pd.read_csv(DATA / "donors.tsv", sep="\t", index_col=0, keep_default_na=False)
D = D[D.modeling_eligible.astype(str) == "True"].copy()
D["stratum"] = D.dx + "|" + D.cohort
sgkf = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=SEED)
tr, te = next(sgkf.split(D, D.stratum, groups=D.split_group))
D["partition"] = "dev"; D.iloc[te, D.columns.get_loc("partition")] = "test"
assert D.groupby("split_group").partition.nunique().max() == 1
D[["dx", "cohort", "batch", "split_group", "affected", "partition", "gsms"]].to_csv(DATA / "split_locked.tsv", sep="\t")
tab = D.groupby(["partition", "cohort", "dx"]).size().unstack(fill_value=0)
print(tab); print(D.groupby("partition").size(), "\ntest fraction per stratum:\n", D.groupby("stratum").partition.apply(lambda v: round((v == "test").mean(), 2)).to_string())
out = {"versions": versions(), "sklearn_version": sklearn.__version__, "seed": SEED,
       "method": "StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=SEED); first fold = test; y = dx|cohort; groups = split_group",
       "n_donors": {k: int(v) for k, v in D.partition.value_counts().items()},
       "n_samples_arrays": {p: int(D[D.partition == p].gsms.str.count(";").add(1).sum()) for p in ("dev", "test")},
       "donors_by_partition_cohort_dx": {p: {c: {k: int(v) for k, v in r.items()} for c, r in tab.loc[p].iterrows()} for p in ("dev", "test")},
       "donors_by_partition_dx": {p: {k: int(v) for k, v in D[D.partition == p].dx.value_counts().items()} for p in ("dev", "test")},
       "affected_by_partition": {p: int((D[D.partition == p].affected.astype(str) == "True").sum()) for p in ("dev", "test")},
       "test_donor_ids": sorted(D[D.partition == "test"].index.tolist())}
write_json(out, RES / "test_split.json")
