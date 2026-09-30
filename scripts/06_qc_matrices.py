"""QC-only gene-level matrices from deposited series matrices (NOT the modeling features).
Used for label-free QC: sex check and donor identity. Output: data/qc/<GSE>[_<GPL>].pkl"""
import glob, pathlib
import numpy as np, pandas as pd
from common import DATA, read_series_matrix, to_gene_level

out = DATA / "qc"; out.mkdir(exist_ok=True)
s = pd.read_csv(DATA / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
for f in sorted(glob.glob(str(DATA / "matrix" / "*_series_matrix.txt.gz"))):
    name = pathlib.Path(f).name.replace("_series_matrix.txt.gz", "")
    gse = name.split("-")[0]
    m = read_series_matrix(f)
    gpl = s[s.gsm.isin(m.columns)].platform.iloc[0]
    raw_max = np.nanmax(m.values)
    if raw_max > 100:  # not log-scale
        m = np.log2(m.clip(lower=2 ** -8))  # GSE12649 deposit is per-gene median-normalised (median 1.0): keep ratios <1
    g = to_gene_level(m, gpl)
    g.to_pickle(out / f"{name}.pkl")
    print(name, gpl, "probes", m.shape[0], "genes", g.shape[0], "samples", g.shape[1], "raw max %.1f" % raw_max,
          "median %.2f" % np.nanmedian(g.values))
