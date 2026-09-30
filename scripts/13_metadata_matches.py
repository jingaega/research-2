"""Metadata-based cross-series donor matching (no expression). For every pair of series that report
sex, age, PMI and pH, count donor pairs with identical sex and age and |dPMI|<=0.5 h and |dpH|<=0.02.
Chance coincidences are estimated by repeating the count after shuffling the age column of series B
(100 permutations). Also compares GSE17612 vs GSE21935 donor ID strings.
Output: results/metadata_matches.json, results/metadata_matches.tsv"""
import itertools
import numpy as np, pandas as pd
from common import DATA, RES, write_json, versions

s = pd.read_csv(DATA / "samples.tsv", sep="\t", keep_default_na=False, na_values=[""], dtype={"donor_local": str})
s = s[s.sample_kind == "tissue"]
cols = ["sex_reported", "age", "pmi", "ph", "dx"]
d = s.dropna(subset=["age", "pmi", "ph"]).groupby(["gse", "donor_local"])[cols].first().reset_index()
d = d[d.sex_reported.isin(["M", "F"])]
rng = np.random.default_rng(0)

def count(a, b):
    m = a.merge(b, on=["sex_reported", "age"], suffixes=("_a", "_b"))
    m = m[(abs(m.pmi_a - m.pmi_b) <= 0.5) & (abs(m.ph_a - m.ph_b) <= 0.02)]
    return m

rows, summ = [], {}
gses = sorted(d.gse.unique())
for ga, gb in itertools.combinations(gses, 2):
    a, b = d[d.gse == ga], d[d.gse == gb]
    m = count(a, b)
    null = []
    for _ in range(100):
        bb = b.copy(); bb["age"] = rng.permutation(bb.age.values); null.append(len(count(a, bb)))
    summ[f"{ga}~{gb}"] = {"n_a": len(a), "n_b": len(b), "matches": len(m), "null_mean": float(np.mean(null)),
                          "null_max": int(np.max(null)), "dx_agree": int((m.dx_a == m.dx_b).sum())}
    for _, r in m.iterrows():
        rows.append({"gse_a": ga, "donor_a": r.donor_local_a, "gse_b": gb, "donor_b": r.donor_local_b, "sex": r.sex_reported,
                     "age": r.age, "pmi_a": r.pmi_a, "pmi_b": r.pmi_b, "ph_a": r.ph_a, "ph_b": r.ph_b, "dx_a": r.dx_a, "dx_b": r.dx_b})
M = pd.DataFrame(rows); M.to_csv(RES / "metadata_matches.tsv", sep="\t", index=False)
# GSE17612 (BA10) vs GSE21935 (BA22): same donor-ID scheme (e.g. C002, S014)
a = s[s.gse == "GSE17612"].groupby("donor_local")[["sex_reported", "age", "dx"]].first()
b = s[s.gse == "GSE21935"].groupby("donor_local")[["sex_reported", "age", "dx"]].first()
shared = a.index.intersection(b.index)
agree = [(i, bool((a.loc[i, "sex_reported"] == b.loc[i, "sex_reported"]) and (a.loc[i, "age"] == b.loc[i, "age"]))) for i in shared]
summ["GSE17612~GSE21935_id_strings"] = {"shared_ids": len(shared), "sex_and_age_agree": sum(x for _, x in agree),
                                        "n_GSE17612": len(a), "n_GSE21935": len(b)}
write_json({"versions": versions(), "rule": "identical sex & age, |dPMI|<=0.5h, |dpH|<=0.02; null = 100 age permutations of series B",
            "pairs": summ}, RES / "metadata_matches.json")
for k, v in summ.items():
    if v.get("matches", v.get("shared_ids", 0)) > 0: print(k, v)
