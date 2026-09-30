"""fRMA probeset matrices -> gene-level (Entrez) matrices, one per (series, chip type).
Probesets mapping to >1 gene were dropped when building the maps; multiple probesets per gene are
averaged (a fixed, data-independent rule). For U133 Plus 2 arrays inside GSE92538 (mixed with U133A
arrays), only probesets also present on U133A are used so that both chip types summarise the same probes.
Output: data/expr/<GSE>__<chip>.pkl"""
import glob, pathlib
import numpy as np, pandas as pd
from common import DATA, to_gene_level

CHIP2GPL = {"HG-U133A": "GPL96", "HG-U133_Plus_2": "GPL570", "HuGene-1_0-st-v1": "GPL6244"}
out = DATA / "expr"; out.mkdir(exist_ok=True)
u133a = set(pd.read_csv(DATA / "annot" / "GPL96_map.tsv", sep="\t", dtype=str).PROBEID)
for f in sorted(glob.glob(str(DATA / "frma" / "*.tsv.gz"))):
    name = pathlib.Path(f).name.replace(".tsv.gz", "")
    of = out / f"{name}.pkl"
    if of.exists(): continue
    gse, chip = name.split("__")
    m = pd.read_csv(f, sep="\t", index_col=0); m.index = m.index.astype(str)
    if gse == "GSE92538" and chip == "HG-U133_Plus_2":
        m = m.loc[m.index.intersection(u133a)]
    gpl = CHIP2GPL.get(chip)
    if gpl is None:
        gpl = "GPL6244" if "hugene" in chip.lower() else None
    g = to_gene_level(m, gpl)
    g.to_pickle(of)
    print(name, gpl, "probesets", m.shape[0], "genes", g.shape[0], "arrays", g.shape[1], "median %.2f" % np.nanmedian(g.values))
