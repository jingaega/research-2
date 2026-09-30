"""Shared helpers: versions, paths, gene-level matrix loading."""
import gzip, io, json, pathlib, platform, sys
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"

def versions():
    import sklearn, scipy
    v = {"python": sys.version.split()[0], "platform": platform.platform(),
         "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
         "scikit_learn": sklearn.__version__}
    try:
        import subprocess
        out = subprocess.run(["Rscript", "-e", 'cat(R.version.string, as.character(packageVersion("limma")), '
                              'as.character(packageVersion("frma")), as.character(packageVersion("affy")), sep="|")'],
                             capture_output=True, text=True, timeout=60).stdout.split("|")
        v.update({"R": out[0], "limma": out[1], "frma": out[2], "affy": out[3]})
    except Exception as e:  # R absent
        v["R"] = f"unavailable: {e}"
    return v

def write_json(obj, path):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp"); tmp.write_text(json.dumps(obj, indent=1, default=_default)); tmp.replace(path)

def _default(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (np.bool_,)): return bool(o)
    return str(o)

def probe_map(gpl):
    if gpl in ("GPL10526", "GPL17027"):  # Brainarray ENTREZG: "<entrez>_at"
        s = pd.read_csv(DATA / "annot" / "entrez_symbol.tsv", sep="\t", dtype=str)
        return None, s
    m = pd.read_csv(DATA / "annot" / f"{gpl}_map.tsv", sep="\t", dtype=str)
    return m.set_index("PROBEID")["ENTREZID"], None

def read_series_matrix(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith("!")]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")

def to_gene_level(df, gpl):
    """Probe -> Entrez gene; multi-probe genes collapsed by mean (data-independent rule)."""
    pm, _ = probe_map(gpl)
    if pm is None:
        ent = df.index.str.replace("_at", "", regex=False)
        ok = np.asarray(pd.Series(ent).str.fullmatch(r"\d+"))
        g = df[ok].copy(); g.index = ent[ok]
    else:
        g = df.loc[df.index.intersection(pm.index)].copy(); g.index = pm.loc[g.index].values
    return g.groupby(level=0).mean()

def symbols():
    s = pd.read_csv(DATA / "annot" / "entrez_symbol.tsv", sep="\t", dtype=str)
    return s.set_index("ENTREZID")
