"""Expression-based sex check (label-free). Y-linked genes vs XIST, per series matrix.
Male call: mean z(Y genes) - z(XIST) > 0 after centring on the midpoint between the two 1-D k-means modes.
A sample is 'ambiguous' when Y-score and XIST disagree (both high or both low)."""
import glob, pathlib, json
import numpy as np, pandas as pd
from common import DATA, RES, symbols, write_json, versions

Y = ["RPS4Y1", "DDX3Y", "KDM5D", "USP9Y", "EIF1AY", "UTY", "ZFY", "NLGN4Y"]
sym = symbols()["SYMBOL"]
s = pd.read_csv(DATA / "samples.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("gsm")

def kmeans1d(x, it=50):
    c = np.percentile(x, [10, 90]).astype(float)
    for _ in range(it):
        lab = np.abs(x[:, None] - c[None]).argmin(1)
        c = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in (0, 1)])
    return c

rows = []
for f in sorted(glob.glob(str(DATA / "qc" / "*.pkl"))):
    g = pd.read_pickle(f); g.index = [sym.get(i, i) for i in g.index]
    g = g.groupby(level=0).mean()
    yg = [x for x in Y if x in g.index]
    reg = s.reindex(g.columns)["region"].fillna("").values
    for rg in sorted(set(reg)):
      gg = g.loc[:, reg == rg]
      if gg.shape[1] < 4: continue
      z = gg.sub(gg.mean(axis=1), axis=0).div(gg.std(axis=1) + 1e-9, axis=0)
      ys = z.loc[yg].mean()
      xs = z.loc["XIST"] if "XIST" in z.index else pd.Series(np.nan, index=z.columns)
      score = (ys - xs.fillna(0)).values
      c = kmeans1d(score); mid = c.mean(); sep = abs(c[1] - c[0])
      for gsm, sc, yv, xv in zip(gg.columns, score, ys.values, xs.values):
        call = "M" if sc > mid else "F"
        # ambiguous: Y-score and XIST agree in sign (both high: XXY/contamination; both low), or score near midpoint
        ambiguous = bool((not np.isnan(xv) and (yv > 0) == (xv > 0)) or abs(sc - mid) < 0.15 * sep)
        rep = s.loc[gsm, "sex_reported"] if gsm in s.index else ""
        rows.append(dict(matrix=pathlib.Path(f).stem, qc_region=rg, gsm=gsm, y_genes=len(yg), y_score=yv, xist_z=xv,
                         score=sc - mid, mode_separation=sep, sex_expr=call, ambiguous=ambiguous, sex_reported=rep,
                         mismatch=bool(rep in ("M", "F") and rep != call)))
r = pd.DataFrame(rows).merge(s[["gse", "title", "dx", "region", "donor_local", "sample_kind"]], left_on="gsm", right_index=True)
r.to_csv(RES / "sex_check.tsv", sep="\t", index=False)
t = r[r.sample_kind == "tissue"]
print(t.groupby(["matrix", "qc_region"]).agg(n=("gsm", "size"), reported=("sex_reported", lambda v: (v != "").sum()),
      mismatches=("mismatch", "sum"), ambiguous=("ambiguous", "sum"), sep=("mode_separation", "first")).to_string())
print(t[t.mismatch | t.ambiguous][["gse", "gsm", "title", "region", "dx", "donor_local", "sex_reported", "sex_expr",
      "y_score", "xist_z", "ambiguous"]].to_string())
