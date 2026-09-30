"""Assemble donor x gene matrices from fRMA gene-level data (no fitting of any kind).
Genes: Entrez genes present on all modeling platforms, minus chrY genes and XIST/TSIX (sex-determining
genes excluded so classifiers cannot use donor sex composition). Replicate arrays of a donor are averaged.
DEV and TEST are written to separate files; the test file is not read until the single final evaluation.
Output: data/features_dev.pkl, data/features_test.pkl (DataFrame donors x genes), data/features_meta.json"""
import glob, pathlib, json
import numpy as np, pandas as pd
from common import DATA, symbols, write_json

sym = symbols()
split = pd.read_csv(DATA / "split_locked.tsv", sep="\t", index_col=0, keep_default_na=False)
arrays = {}
for f in sorted(glob.glob(str(DATA / "expr" / "*.pkl"))):
    g = pd.read_pickle(f)
    for c in g.columns: arrays[c] = g[c]
needed = sorted({x for v in split.gsms for x in v.split(";")})
series_genes = [set(pd.read_pickle(f).index) for f in glob.glob(str(DATA / "expr" / "*.pkl"))]
common = set.intersection(*series_genes)
chrom = sym.reindex(sorted(common))["CHR"].astype(str)
drop = set(chrom[chrom == "Y"].index) | set(sym[sym.SYMBOL.isin(["XIST", "TSIX"])].index)
genes = sorted(common - drop)
rows = {}
for did, r in split.iterrows():
    rows[did] = pd.concat([arrays[x].reindex(genes) for x in r.gsms.split(";")], axis=1).mean(axis=1)
X = pd.DataFrame(rows).T
assert X.notna().all().all()
X.loc[split.index[split.partition == "dev"]].to_pickle(DATA / "features_dev.pkl")
X.loc[split.index[split.partition == "test"]].to_pickle(DATA / "features_test.pkl")
write_json({"n_genes": len(genes), "n_common_before_sex_filter": len(common), "n_sex_genes_dropped": len(drop & common),
            "n_dev": int((split.partition == "dev").sum()), "n_test": int((split.partition == "test").sum())}, DATA / "features_meta.json")
print("genes:", len(genes), "(dropped sex genes:", len(drop & common), ") dev:", (split.partition == "dev").sum(), "test:", (split.partition == "test").sum())
