"""Sex check on fRMA gene-level data (Toker et al. 2016 style): separate two-cluster calls for the
Y-gene score and for XIST, computed per (series, chip, region). Label-free.
 - sex_expr: 'M' if Y-call male and XIST-call low; 'F' if Y-call female and XIST-call high;
   'AMBIG' if the two calls disagree (e.g. XXY, contamination/mixture, or failed probes).
 - HuGene (GSE35978): XIST transcript cluster is multi-mapped and excluded from the map, so the call
   uses the Y score only (XIST taken from the raw transcript cluster if present).
Output: results/sex_check.tsv (+ summary in DATA_AUDIT.json later)."""
import glob, pathlib
import numpy as np, pandas as pd
from common import DATA, RES, symbols

Y = ["RPS4Y1", "DDX3Y", "KDM5D", "USP9Y", "EIF1AY", "UTY", "ZFY", "NLGN4Y"]
sym = symbols()["SYMBOL"]
s = pd.read_csv(DATA / "samples.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("gsm")

def two_means(x, it=100):
    x = np.asarray(x, float); c = np.percentile(x, [5, 95])
    for _ in range(it):
        lab = np.abs(x[:, None] - c[None]).argmin(1)
        c = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in (0, 1)])
    lab = np.abs(x[:, None] - c[None]).argmin(1)
    return lab, c

def call(x, min_sep):
    """1 = high cluster. If the two cluster centres are closer than min_sep (single-sex group),
    fall back to an absolute threshold taken from the other datasets (returned as None)."""
    lab, c = two_means(x)
    return (lab == c.argmax()).astype(int), abs(c[1] - c[0])

# Collect absolute Y-score / XIST per sample, then call per CHIP TYPE pooled over all series and regions
# (per-unit clustering is wrong for single-sex series such as GSE54567/GSE54572 = males only).
recs = []
for f in sorted(glob.glob(str(DATA / "expr" / "*.pkl"))):
    name = pathlib.Path(f).stem; gse, chip = name.split("__")
    g = pd.read_pickle(f)
    g.index = [sym.get(i, i) for i in g.index]; g = g.groupby(level=0).mean()
    yg = [x for x in Y if x in g.index]
    meta = s.reindex(g.columns)
    cols = meta.index[meta.sample_kind == "tissue"]
    ysc = g.loc[yg, cols].mean()
    xs = g.loc["XIST", cols] if "XIST" in g.index else pd.Series(np.nan, index=cols)
    for gsm in cols:
        recs.append(dict(expr=name, chip=chip, gse=gse, gsm=gsm, y_score=float(ysc[gsm]),
                         xist=float(xs[gsm]) if not np.isnan(xs[gsm]) else np.nan))
R0 = pd.DataFrame(recs)
rows = []
for chip, d in R0.groupby("chip"):
    ycall, yc = two_means(d.y_score.values); ymale = ycall == yc.argmax(); ysep = abs(yc[1] - yc[0])
    if d.xist.notna().all():
        xcall, xc = two_means(d.xist.values); xfem = xcall == xc.argmax(); xsep = abs(xc[1] - xc[0])
    else:
        xfem = np.full(len(d), np.nan); xsep = np.nan
    ymid = yc.mean()
    for i, (_, r) in enumerate(d.iterrows()):
        m = s.loc[r.gsm]
        if np.isnan(xsep): sx = "M" if ymale[i] else "F"
        elif ymale[i] and not xfem[i]: sx = "M"
        elif (not ymale[i]) and xfem[i]: sx = "F"
        else: sx = "AMBIG"
        rep = m.sex_reported
        rows.append(dict(expr=r.expr, chip=chip, gse=r.gse, region=m.region, gsm=r.gsm, title=m.title, dx=m.dx,
                         donor_local=m.donor_local, y_score=round(r.y_score, 3),
                         xist=None if np.isnan(r.xist) else round(r.xist, 3),
                         y_margin=round((r.y_score - ymid) / ysep, 3), y_sep=round(ysep, 2),
                         xist_sep=None if np.isnan(xsep) else round(xsep, 2), sex_reported=rep, sex_expr=sx,
                         mismatch=bool(rep in ("M", "F") and sx in ("M", "F") and rep != sx)))
r = pd.DataFrame(rows)
r.to_csv(RES / "sex_check.tsv", sep="\t", index=False)
print(r.groupby(["expr", "region"]).agg(n=("gsm", "size"), y_sep=("y_sep", "first"), xist_sep=("xist_sep", "first"),
      mismatch=("mismatch", "sum"), ambiguous=("sex_expr", lambda v: (v == "AMBIG").sum()),
      unreported=("sex_reported", lambda v: (~v.isin(["M", "F"])).sum())).to_string())
bad = r[r.mismatch | (r.sex_expr == "AMBIG")]
print(bad[["gse", "region", "gsm", "title", "dx", "donor_local", "sex_reported", "sex_expr", "y_score", "xist"]].to_string())
