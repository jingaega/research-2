"""Expression-based donor identity audit (label-free; diagnosis is only used afterwards to REPORT agreement).

1. Units = (expression file, region). Within each unit genes are z-scored across its tissue samples.
2. Identity genes: genes whose per-donor deviations are reproducible across regions/chips of the SAME donor
   (genotype-driven expression, e.g. deletion polymorphisms / strong cis-eQTLs). Selected on calibration
   pairs with known identity from metadata: Sibille BA9<->BA25 pair numbers, GSE53987 metadata donors across
   PFC/HPC/STR, GSE92538 subject ids across replicate arrays. Sex-chromosome genes excluded.
3. Validation on held-out known pairs: GSE35978 cerebellum <-> parietal (metadata donor keys).
4. Scan every pair of units (and each unit against itself for hidden duplicates): Pearson correlation on
   identity genes; a candidate match is a MUTUAL best hit whose robust z-score against all correlations of that
   unit pair is >= Z_MATCH. Confidence tiers are defined from the validation set.
Outputs: results/identity_genes.tsv, results/identity_validation.json, results/identity_matches.tsv"""
import glob, itertools, json, pathlib
import numpy as np, pandas as pd
from common import DATA, RES, symbols, write_json, versions

sym = symbols()
s = pd.read_csv(DATA / "samples.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("gsm")
units = {}
for f in sorted(glob.glob(str(DATA / "expr" / "*.pkl"))):
    name = pathlib.Path(f).stem; g = pd.read_pickle(f)
    meta = s.reindex(g.columns)
    for rg in sorted(meta.region.dropna().unique()):
        cols = meta.index[(meta.region == rg) & (meta.sample_kind == "tissue")]
        if len(cols) < 8: continue
        units[f"{name}|{rg}"] = g[cols]
common = sorted(set.intersection(*[set(u.index) for u in units.values()]))
chrom = sym.reindex(common)["CHR"].astype(str)
auto = [g for g in common if chrom.get(g, "") not in ("X", "Y", "nan", "")]
Z = {k: (u.loc[auto].sub(u.loc[auto].mean(axis=1), axis=0).div(u.loc[auto].std(axis=1) + 1e-9, axis=0)) for k, u in units.items()}
print("units:", len(units), "| genes common to all units:", len(common), "| autosomal:", len(auto))

def donor_key(gsm):
    r = s.loc[gsm]; return f"{r.gse}:{r.donor_local}"

def known_pairs(unit_names):
    """All (unitA, gsmA, unitB, gsmB) with the same metadata donor key in different units."""
    idx = [(u, c, donor_key(c)) for u in unit_names for c in Z[u].columns]
    by = {}
    for u, c, k in idx: by.setdefault(k, []).append((u, c))
    pairs = []
    for k, lst in by.items():
        for (u1, c1), (u2, c2) in itertools.combinations(lst, 2):
            if u1 != u2: pairs.append((u1, c1, u2, c2))
    return pairs

def sibille_pairs():
    out = []
    for a, b in (("GSE54567", "GSE54572"), ("GSE54568", "GSE54571")):
        ua = [u for u in Z if u.startswith(a)][0]; ub = [u for u in Z if u.startswith(b)][0]
        mb = {s.loc[c, "donor_local"]: c for c in Z[ub].columns}
        out += [(ua, c, ub, mb[s.loc[c, "donor_local"]]) for c in Z[ua].columns if s.loc[c, "donor_local"] in mb]
    return out

def pairs_for(prefix):
    return known_pairs([u for u in Z if u.startswith(prefix)])

def excl_ambiguous(pairs, gse):
    # GSE35978 metadata keys that hold >1 sample of the same region are not trustworthy identities
    bad = set()
    for u in [u for u in Z if u.startswith(gse)]:
        v = pd.Series([donor_key(c) for c in Z[u].columns]).value_counts(); bad |= set(v[v > 1].index)
    return [p for p in pairs if donor_key(p[1]) not in bad]

calib = {"sibille_BA9_vs_BA25": sibille_pairs(),
         "GSE53987_regions": pairs_for("GSE53987"),
         "GSE92538_replicate_arrays": pairs_for("GSE92538")}
valid = {"GSE35978_CB_vs_PARIETAL": excl_ambiguous(pairs_for("GSE35978"), "GSE35978")}
print({k: len(v) for k, v in {**calib, **valid}.items()})

def gene_r(pairs):
    A = np.stack([Z[u1][c1].values for u1, c1, _, _ in pairs], 1); B = np.stack([Z[u2][c2].values for _, _, u2, c2 in pairs], 1)
    A = (A - A.mean(1, keepdims=True)) / (A.std(1, keepdims=True) + 1e-9); B = (B - B.mean(1, keepdims=True)) / (B.std(1, keepdims=True) + 1e-9)
    return (A * B).mean(1)

R = pd.DataFrame({k: gene_r(v) for k, v in {**calib, **valid}.items()}, index=auto)
R["calib_min"] = R[list(calib)].min(axis=1); R["calib_mean"] = R[list(calib)].mean(axis=1)
N_ID = 300
idg = R.sort_values("calib_min", ascending=False).index[:N_ID].tolist()
R.assign(symbol=sym.reindex(R.index)["SYMBOL"], selected=R.index.isin(idg)).sort_values("calib_min", ascending=False) \
 .to_csv(RES / "identity_genes.tsv", sep="\t")
print("top identity genes:", sym.reindex(idg[:25])["SYMBOL"].tolist())
print("validation-set gene r: selected mean %.3f vs all-genes mean %.3f" % (R.loc[idg, "GSE35978_CB_vs_PARIETAL"].mean(), R["GSE35978_CB_vs_PARIETAL"].mean()))

def scan(ua, ub, genes=idg):
    A = Z[ua].loc[genes]; B = Z[ub].loc[genes]
    C = np.corrcoef(A.T.values, B.T.values)[:A.shape[1], A.shape[1]:]
    same = ua == ub
    if same: np.fill_diagonal(C, np.nan)
    flat = C[~np.isnan(C)]; med = np.median(flat); mad = 1.4826 * np.median(np.abs(flat - med)) + 1e-9
    best_b = np.nanargmax(C, 1); best_a = np.nanargmax(C, 0)
    out = []
    for i in range(C.shape[0]):
        j = best_b[i]
        if best_a[j] != i: continue
        if same and j < i: continue
        row = np.sort(C[i][~np.isnan(C[i])])[::-1]
        out.append(dict(unit_a=ua, gsm_a=A.columns[i], unit_b=ub, gsm_b=B.columns[j], r=C[i, j],
                        z=(C[i, j] - med) / mad, r_second=row[1] if len(row) > 1 else np.nan))
    return pd.DataFrame(out)

# validation: how well do known CB<->PARIETAL pairs get recovered?
va = valid["GSE35978_CB_vs_PARIETAL"]
ucb = [u for u in Z if u.startswith("GSE35978") and u.endswith("|CB")][0]; upa = [u for u in Z if u.startswith("GSE35978") and u.endswith("|PARIETAL")][0]
sc = scan(ucb, upa); truth = {(a, b) for _, a, _, b in [(p[0], p[1], p[2], p[3]) if p[0] == ucb else (p[2], p[3], p[0], p[1]) for p in va]}
sc["true"] = [(a, b) in truth for a, b in zip(sc.gsm_a, sc.gsm_b)]
vres = {"n_known_pairs": len(truth), "n_mutual_best": int(len(sc))}
for zt in (3, 4, 5, 6, 8):
    called = sc[sc.z >= zt]
    vres[f"z>={zt}"] = {"called": int(len(called)), "true_positive": int(called.true.sum()),
                        "false_positive": int((~called.true).sum()), "sensitivity": round(called.true.sum() / len(truth), 3)}
print("validation:", json.dumps(vres))
Z_MATCH = 5.0
allm = []
names = sorted(Z)
for ua, ub in itertools.combinations_with_replacement(names, 2):
    m = scan(ua, ub)
    if len(m): allm.append(m[m.z >= 3])
M = pd.concat(allm, ignore_index=True)
for side in ("a", "b"):
    M[f"gse_{side}"] = s.reindex(M[f"gsm_{side}"])["gse"].values
    for c in ("dx", "region", "sex_reported", "age", "pmi", "ph", "donor_local", "title"):
        M[f"{c}_{side}"] = s.reindex(M[f"gsm_{side}"])[c].values
M["same_metadata_donor"] = (M.gse_a == M.gse_b) & (M.donor_local_a == M.donor_local_b)
M["dx_agree"] = M.dx_a == M.dx_b
M.to_csv(RES / "identity_matches.tsv", sep="\t", index=False)
write_json({"versions": versions(), "n_identity_genes": N_ID, "n_common_genes": len(common), "units": {k: int(v.shape[1]) for k, v in Z.items()},
            "calibration_pairs": {k: len(v) for k, v in calib.items()}, "validation": vres, "z_match": Z_MATCH},
           RES / "identity_validation.json")
new = M[(M.z >= Z_MATCH) & ~M.same_metadata_donor]
print("\nmatches z>=%.0f not explained by within-series metadata donor keys, by series pair:" % Z_MATCH)
print(new.groupby(["gse_a", "gse_b"]).agg(n=("r", "size"), dx_agree=("dx_agree", "sum"), median_z=("z", "median")).to_string())
