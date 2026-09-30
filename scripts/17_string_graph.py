"""Curated knowledge graph: STRING v12 human, combined_score >= 700 (high confidence), mapped to Entrez via
STRING preferred_name -> org.Hs.eg.db symbol. Restricted to modeling genes. Undirected, unweighted.
Output: data/graph_string.npz (edge index over feature gene order)"""
import gzip
import numpy as np, pandas as pd
from common import DATA, symbols, write_json

genes = pd.read_pickle(DATA / "features_dev.pkl").columns.tolist()
gi = {g: i for i, g in enumerate(genes)}
sym = symbols().reset_index()
s2e = sym.dropna(subset=["SYMBOL"]).drop_duplicates("SYMBOL").set_index("SYMBOL")["ENTREZID"]
info = pd.read_csv(DATA / "knowledge" / "9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1])
info.columns = ["protein", "name"]; p2e = info.set_index("protein")["name"].map(s2e).dropna()
L = pd.read_csv(DATA / "knowledge" / "9606.protein.links.v12.0.txt.gz", sep=" ")
L = L[L.combined_score >= 700]
a = L.protein1.map(p2e); b = L.protein2.map(p2e)
ok = a.notna() & b.notna()
e = pd.DataFrame({"a": a[ok].map(gi), "b": b[ok].map(gi)}).dropna().astype(int)
e = e[e.a != e.b]
e = pd.DataFrame(np.sort(e.values, axis=1), columns=["a", "b"]).drop_duplicates()
np.savez_compressed(DATA / "graph_string.npz", edges=e.values, n=len(genes))
deg = np.bincount(e.values.ravel(), minlength=len(genes))
write_json({"source": "STRING v12.0 9606 protein.links, combined_score>=700", "n_nodes": len(genes), "n_edges": int(len(e)),
            "n_isolated_nodes": int((deg == 0).sum()), "median_degree": float(np.median(deg)), "mean_degree": float(deg.mean())},
           DATA / "graph_string_meta.json")
print("STRING edges:", len(e), "isolated:", (deg == 0).sum(), "median deg:", np.median(deg), "mean deg: %.1f" % deg.mean())
