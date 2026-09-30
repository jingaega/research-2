"""Build the donor-deduplicated table, modeling-sample selection, class inclusion and DATA_AUDIT.json.
No expression modeling. Diagnosis labels are used only as metadata (counts, agreement checks)."""
import json, collections
import numpy as np, pandas as pd
from common import DATA, RES, ROOT, write_json, versions

s = pd.read_csv(DATA / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
sex = pd.read_csv(RES / "sex_check.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("gsm")
M = pd.read_csv(RES / "identity_matches_pc10.tsv", sep="\t", dtype=str, keep_default_na=False)
M["z"] = M.z.astype(float)
idv = json.load(open(RES / "identity_validation_pc10.json"))
meta_m = json.load(open(RES / "metadata_matches.json"))
Z_MATCH, Z_PROB, Z_SPLIT = idv["z_match"], idv["z_probable"], 3.0

MODEL_SERIES = ["GSE12649", "GSE35978", "GSE17612", "GSE21138", "GSE53987", "GSE54567", "GSE54568", "GSE54571", "GSE54572", "GSE92538"]
BANK = {"GSE12649": "Stanley", "GSE35978": "Stanley", "GSE17612": "CharingCross", "GSE21138": "Victoria",
        "GSE53987": "Pittsburgh", "GSE54567": "Pittsburgh", "GSE54568": "Pittsburgh", "GSE54571": "Pittsburgh",
        "GSE54572": "Pittsburgh", "GSE92538": "Pritzker"}
KEEP_DX = {"SCZ", "BD", "MDD", "CTL"}

# ---------------- sample-level exclusions (before donor building) ----------------
n = s[s.gse.isin(MODEL_SERIES) & (s.sample_kind == "tissue")].copy()
n["excluded_reason"] = ""
n.loc[~n.dx.isin(KEEP_DX), "excluded_reason"] = "diagnosis not SCZ/BD/MDD/CTL (e.g. 'bipolar (not bipolar)')"
n.loc[(n.gse == "GSE53987") & (n.region == "HPC"), "excluded_reason"] = \
    "GSE53987 hippocampus: expression identity shows HPC arrays belong to other recorded donors (label mix-up)"
n["sex_expr"] = sex.reindex(n.gsm)["sex_expr"].values
n["bank"] = n.gse.map(BANK)
n["node"] = n.gsm
use = n[n.excluded_reason == ""].copy()

# ---------------- union-find over samples ----------------
par = {g: g for g in use.gsm}
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
def union(a, b):
    if a in par and b in par: par[find(a)] = find(b)

# (1) within-series metadata donor keys (GSE92538 subject id, Sibille pair label, GSE12649/17612/21138 IDs,
#     GSE35978/GSE53987 identical (dx,sex,age,PMI,pH[,race]) tuples)
for (g, d), grp in use.groupby(["gse", "donor_local"]):
    if d == "": continue
    x = grp.gsm.tolist()
    for y in x[1:]: union(x[0], y)
# Sibille: same pair label across BA9/BA25 series of the same sex = same donor (verified by expression: 42/42 mutual-best
# hits between same-sex BA9/BA25 series carry identical pair labels)
sib = use[use.gse.isin(["GSE54567", "GSE54568", "GSE54571", "GSE54572"])]
for d, grp in sib.groupby("donor_local"):
    x = grp.gsm.tolist()
    for y in x[1:]: union(x[0], y)
meta_union = {g: find(g) for g in use.gsm}

# (2) expression identity edges
E = M[M.gsm_a.isin(use.gsm) & M.gsm_b.isin(use.gsm) & (M.bank_a == M.bank_b)].copy()
E["sex_a"] = sex.reindex(E.gsm_a)["sex_expr"].values; E["sex_b"] = sex.reindex(E.gsm_b)["sex_expr"].values
E["sex_agree"] = E.sex_a == E.sex_b
E["dx_agree"] = E.dx_a == E.dx_b
E["confirmed"] = (E.z >= Z_MATCH) | ((E.z >= Z_PROB) & E.dx_agree & E.sex_agree)
E["split_edge"] = E.confirmed | ((E.z >= Z_SPLIT) & E.dx_agree & E.sex_agree)
conf = E[E.confirmed]
for a, b in zip(conf.gsm_a, conf.gsm_b): union(a, b)
donor_of = {g: find(g) for g in use.gsm}
# split groups: confirmed + liberal edges
spar = dict(donor_of)
def sfind(x):
    while spar[x] != x: spar[x] = spar[spar[x]]; x = spar[x]
    return x
for g in list(spar): spar[g] = donor_of[g]
roots = {r: r for r in set(donor_of.values())}
def rfind(x):
    while roots[x] != x: roots[x] = roots[roots[x]]; x = roots[x]
    return x
for a, b in zip(E[E.split_edge].gsm_a, E[E.split_edge].gsm_b):
    ra, rb = rfind(donor_of[a]), rfind(donor_of[b])
    if ra != rb: roots[ra] = rb
use["donor_uf"] = use.gsm.map(donor_of)
use["split_uf"] = use.donor_uf.map(rfind)
# stable IDs
dmap = {r: f"D{i:04d}" for i, r in enumerate(sorted(use.donor_uf.unique()))}
gmap = {r: f"G{i:04d}" for i, r in enumerate(sorted(use.split_uf.unique()))}
use["donor_id"] = use.donor_uf.map(dmap); use["split_group"] = use.split_uf.map(gmap)

# ---------------- donor-level attributes ----------------
def agg(grp):
    dxs = sorted(set(grp.dx)); sr = sorted(set(grp.sex_reported) - {""}); se = sorted(set(grp.sex_expr) - {""})
    se_clear = sorted(set(se) & {"M", "F"})
    return pd.Series({
        "series": ",".join(sorted(set(grp.gse))), "bank": ",".join(sorted(set(grp.bank))), "n_samples": len(grp),
        "regions": ",".join(sorted(set(grp.region))), "dx_labels": ",".join(dxs), "label_conflict": len(dxs) > 1,
        "sex_reported": ",".join(sr), "sex_expr": ",".join(se),
        "sex_mismatch": bool(sr and se_clear and (len(sr) > 1 or set(se_clear) != set(sr))),
        "both_sexes_in_expression": len(se_clear) > 1, "mixed_sample": "MIXED" in se,
        "split_group": grp.split_group.iloc[0]})
D = use.groupby("donor_id").apply(agg, include_groups=False)
D["cohort"] = D.bank
assert (~D.bank.str.contains(",")).all(), "donor spans banks"

# ---------------- modeling sample selection (one cortical sample per donor) ----------------
# Stanley: GSE35978 parietal only (GSE12649 audit-only: its unmatched donors are presumed undetected duplicates)
# Pittsburgh: GSE53987 PFC > Sibille BA9 (BA25 and STR not used); CharingCross: GSE17612 BA10; Victoria: GSE21138 BA46;
# Pritzker: average of the subject's U133A arrays, else its U133 Plus 2 arrays (U133A probes).
PREF = [("GSE35978", "PARIETAL"), ("GSE53987", "PFC_BA46"), ("GSE54567", "PFC_BA9"), ("GSE54568", "PFC_BA9"),
        ("GSE17612", "PFC_BA10"), ("GSE21138", "PFC_BA46"), ("GSE92538", "PFC_DLPFC")]
sel = []
for did, grp in use.groupby("donor_id"):
    for g, rg in PREF:
        cand = grp[(grp.gse == g) & (grp.region == rg)]
        if len(cand):
            if g == "GSE92538":
                a = cand[cand.chip == "U133A"]; cand = a if len(a) else cand
                batch = "GSE92538_U133A" if len(a) else "GSE92538_U133Plus2(U133Aprobes)"
            else:
                batch = g
            sel.append({"donor_id": did, "gsms": ";".join(cand.gsm), "series": g, "region": rg, "batch": batch,
                        "n_arrays_averaged": len(cand)})
            break
S = pd.DataFrame(sel).set_index("donor_id")
D = D.join(S, rsuffix="_model")
# Label = the modeling sample's own recorded diagnosis. A different label on another sample of the same donor
# (audit-only series, other region) marks the donor 'affected'. If the modeling samples themselves carry
# different labels (e.g. two Pritzker subject IDs that are one person), the label is unresolvable -> not eligible.
gsm_dx = use.set_index("gsm").dx
def model_dx(g):
    if not isinstance(g, str): return ""
    labs = set(gsm_dx.reindex(g.split(";")))
    return labs.pop() if len(labs) == 1 else "UNRESOLVABLE"
D["dx"] = D.gsms.map(model_dx)
D["modeling_eligible"] = D.gsms.notna() & D.dx.isin(KEEP_DX)
gsm_sex = use.set_index("gsm").sex_expr
D["model_sample_sex_expr"] = D.gsms.map(lambda g: ",".join(sorted(set(gsm_sex.reindex(g.split(";"))))) if isinstance(g, str) else "")
D["affected"] = D.label_conflict | D.sex_mismatch | D.both_sexes_in_expression | D.mixed_sample
def reason(r):
    out = []
    if r.label_conflict: out.append("identity match to a sample with a different diagnosis (%s)" % r.dx_labels)
    if r.sex_mismatch: out.append("expression sex %s vs reported %s" % (r.sex_expr, r.sex_reported))
    if r.both_sexes_in_expression: out.append("samples of this donor show both sexes")
    if r.mixed_sample: out.append("a sample shows Y and XIST both high (mixture/XXY)")
    return "; ".join(out)
D["affected_reason"] = D.apply(reason, axis=1)

# ---------------- class inclusion (pre-declared: >= 3 cohorts with >= 8 unique donors) ----------------
MIN_DONORS = 8
elig = D[D.modeling_eligible]
cnt = elig[elig.dx.isin(["SCZ", "BD", "MDD"])].groupby(["dx", "cohort"]).size().unstack(fill_value=0)
cnt_clean = elig[elig.dx.isin(["SCZ", "BD", "MDD"]) & ~elig.affected].groupby(["dx", "cohort"]).size().unstack(fill_value=0)
inclusion = {}
for dx in ["SCZ", "BD", "MDD"]:
    row = cnt.loc[dx] if dx in cnt.index else pd.Series(dtype=int)
    rowc = cnt_clean.loc[dx] if dx in cnt_clean.index else pd.Series(dtype=int)
    contributing = sorted(row[row >= MIN_DONORS].index.tolist())
    contributing_clean = sorted(rowc[rowc >= MIN_DONORS].index.tolist())
    inclusion[dx] = {"donors_per_cohort": {k: int(v) for k, v in row.items()},
                     "donors_per_cohort_excluding_affected": {k: int(v) for k, v in rowc.items()},
                     "contributing_cohorts": contributing, "n_contributing": len(contributing),
                     "n_contributing_excluding_affected": len(contributing_clean),
                     "included": len(contributing) >= 3 and len(contributing_clean) >= 3}
CLASSES = [k for k, v in inclusion.items() if v["included"]]
print("class inclusion:", {k: (v["n_contributing"], v["included"]) for k, v in inclusion.items()})

# ---------------- save ----------------
D.to_csv(DATA / "donors.tsv", sep="\t")
use[["gsm", "gse", "region", "dx", "sex_reported", "sex_expr", "donor_local", "donor_id", "split_group", "bank"]].to_csv(DATA / "sample_donor_map.tsv", sep="\t", index=False)
E.to_csv(RES / "identity_edges_used.tsv", sep="\t", index=False)

# ---------------- audit summaries ----------------
def series_table():
    out = {}
    raw = pd.read_csv(DATA / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
    txt = {}
    for g in sorted(raw.gse.unique()):
        f = DATA / "geo_meta" / f"{g}.txt"
        tt = [l.split(" = ", 1)[1] for l in f.read_text(errors="replace").splitlines() if l.startswith("!Series_title")]
        txt[g] = tt[0] if tt else ""
    decision = {
        "GSE11223": "EXCLUDED: colon biopsies (ulcerative colitis), not brain, not psychiatric",
        "GSE28475": "EXCLUDED from modeling: ASD methods study (frozen/fixed, DASL/IVT, reference RNA); ASD has < 3 cohorts",
        "GSE28521": "EXCLUDED from modeling: ASD has < 3 cohorts (GSE28521, GSE28475 only)",
        "GSE29555": "EXCLUDED from modeling: alcohol use disorder, single series",
        "GSE12649": "AUDIT ONLY: Stanley Array collection; 80/102 donors matched to GSE35978; rest presumed undetected duplicates",
        "GSE35978": "INCLUDED (parietal cortex; cerebellum used only for identity/sex audit)",
        "GSE17612": "INCLUDED", "GSE21138": "INCLUDED",
        "GSE53987": "INCLUDED (PFC only; STR for audit; HPC excluded: label mix-ups)",
        "GSE54567": "INCLUDED (BA9 males)", "GSE54568": "INCLUDED (BA9 females)",
        "GSE54571": "AUDIT ONLY (BA25 females; same donors as GSE54568)", "GSE54572": "AUDIT ONLY (BA25 males; same donors as GSE54567)",
        "GSE92538": "INCLUDED (ADDED by us; not in candidate list)",
        "GSE21935": "SCREENED, NOT ADDED (same Charing Cross cohort as GSE17612, BA22)",
        "GSE5388": "SCREENED, NOT ADDED (Stanley Array collection BD/CTL)"}
    for g in sorted(raw.gse.unique()):
        x = raw[raw.gse == g]; t = x[x.sample_kind == "tissue"]
        out[g] = {"title": txt[g], "platforms": sorted(set(x.platform)), "regions": sorted(set(t.region)),
                  "n_samples_total": int(len(x)), "n_tissue_samples": int(len(t)),
                  "n_non_tissue_samples": {k: int(v) for k, v in x[x.sample_kind != "tissue"].sample_kind.value_counts().items()},
                  "n_donors_within_series": int(t[t.donor_local != ""].donor_local.nunique()) if (t.donor_local != "").any() else None,
                  "donor_id_source": {"GSE35978": "no donor id in GEO; keyed on identical (dx,sex,age,PMI,pH)",
                                      "GSE53987": "no donor id in GEO; keyed on identical (dx,sex,age,PMI,pH,race)",
                                      "GSE12649": "number in title (BA46-N)"}.get(g, "title/characteristics"),
                  "diagnosis_counts_samples": {k: int(v) for k, v in t.dx.value_counts().items()},
                  "diagnosis_counts_donors": {k: int(v) for k, v in t[t.donor_local != ""].groupby("dx").donor_local.nunique().items()},
                  "brain_bank": BANK.get(g, {"GSE28521": "Autism Tissue Program (ASD)", "GSE28475": "ASD (methods study)",
                                             "GSE29555": "alcohol study", "GSE11223": "not brain", "GSE21935": "CharingCross",
                                             "GSE5388": "Stanley"}.get(g, "")),
                  "decision": decision.get(g, "")}
    return out

sx = sex.reset_index()
sx = sx[sx.gse.isin(MODEL_SERIES + ["GSE21935", "GSE5388"])]
sex_summary = {g: {"n": int(len(x)), "reported": int(x.sex_reported.isin(["M", "F"]).sum()),
                   "mismatch": int((x.mismatch == "True").sum()), "mixed": int((x.sex_expr == "MIXED").sum()),
                   "low_both": int((x.sex_expr == "LOW_BOTH").sum())} for g, x in sx.groupby("gse")}
sex_list = sx[(sx.mismatch == "True") | sx.sex_expr.isin(["MIXED", "LOW_BOTH"])][
    ["gse", "region", "gsm", "title", "dx", "sex_reported", "sex_expr", "y_score", "xist"]].to_dict("records")
donors_both = D[D.both_sexes_in_expression].reset_index()[["donor_id", "series", "regions", "dx_labels", "sex_reported", "sex_expr"]].to_dict("records")

dup = E[E.confirmed & (E.gse_a != E.gse_b)]
dup_pairs = dup.assign(tier_conf=np.where(dup.z >= Z_MATCH, "high (z >= match threshold, above all negative controls)",
                                          "medium (probable tier + diagnosis and sex agree)"))
cross = collections.Counter()
for _, r in dup_pairs.iterrows(): cross[f"{r.gse_a}~{r.gse_b}"] += 1
multi_series_donors = D[D.series.str.contains(",")]
wide = {"n_donors_total": int(len(D)), "n_donors_in_more_than_one_series": int(len(multi_series_donors)),
        "by_series_combination": {k: int(v) for k, v in multi_series_donors.series.value_counts().items()},
        "by_confidence": {"high": int(multi_series_donors.index.isin(use[use.gsm.isin(dup[dup.z >= Z_MATCH].gsm_a)].donor_id).sum()),
                          "medium_only": None},
        "confirmed_cross_series_sample_pairs": {k: int(v) for k, v in cross.items()},
        "n_split_groups": int(D.split_group.nunique()),
        "n_split_groups_with_more_than_one_donor": int((D.split_group.value_counts() > 1).sum())}
hi_donors = set(use[use.gsm.isin(set(dup[dup.z >= Z_MATCH].gsm_a) | set(dup[dup.z >= Z_MATCH].gsm_b))].donor_id)
wide["by_confidence"] = {"high": int(len(set(multi_series_donors.index) & hi_donors)),
                         "medium_only": int(len(set(multi_series_donors.index) - hi_donors))}

conflicts = E[(E.z >= Z_MATCH) & ~E.dx_agree][["gsm_a", "title_a", "gse_a", "dx_a", "gsm_b", "title_b", "gse_b", "dx_b", "z"]].to_dict("records")
within_unexplained = M[(M.gse_a == M.gse_b) & (M.z >= Z_PROB) & (M.donor_local_a != M.donor_local_b)][
    ["gse_a", "title_a", "region_a", "dx_a", "title_b", "region_b", "dx_b", "z", "tier"]].to_dict("records")

model = D[D.modeling_eligible & D.dx.isin(KEEP_DX)]
tab = lambda d: {c: {k: int(v) for k, v in x.items()} for c, x in d.groupby("cohort").dx.value_counts().unstack(fill_value=0).T.items()} if len(d) else {}
tab = lambda d: d.groupby(["cohort", "dx"]).size().unstack(fill_value=0).astype(int).to_dict(orient="index")

tasks = {}
for dx in CLASSES:
    coh = inclusion[dx]["contributing_cohorts"]
    t = model[model.cohort.isin(coh) & model.dx.isin([dx, "CTL"])]
    tasks[f"A_{dx}_vs_CTL"] = {"cohorts": coh, "classes": [dx, "CTL"], "donors": tab(t),
                               "n_total": int(len(t)), "n_affected": int(t.affected.sum())}
multi = [c for c in model.cohort.unique() if sum(c in inclusion[d]["contributing_cohorts"] for d in CLASSES) >= 2]
tB = model[model.cohort.isin(multi) & model.dx.isin(CLASSES)]
tasks["B_multiclass"] = {"cohorts": sorted(multi), "classes": CLASSES, "donors": tab(tB), "n_total": int(len(tB)),
                         "n_affected": int(tB.affected.sum()),
                         "rule": "cohorts contributing >= 2 included disorders (>= 8 donors each); single-disorder cohorts would make cohort identity = label"}

affected_list = D[D.modeling_eligible & D.affected].reset_index()[["donor_id", "cohort", "series", "dx", "dx_labels", "sex_reported", "sex_expr", "affected_reason"]].to_dict("records")
unresolvable = D[D.dx == "UNRESOLVABLE"].reset_index()[["donor_id", "series", "dx_labels", "gsms"]].to_dict("records")
audit = {
    "created": "2026-09-30", "versions": versions(),
    "scope": "psychiatric postmortem brain microarrays (NCBI GEO); RNA-seq out of scope",
    "candidate_list_status": "the 13 candidate accessions equal the accession list in Gandal et al. 2018 (Science) acknowledgements; GSE11223 is not brain tissue",
    "series": series_table(),
    "preprocessing": {"start": "raw CEL files, frozen RMA (single-array; no cross-sample information)",
                      "gene_mapping": "Bioconductor annotation; multi-gene probesets dropped; probesets averaged per Entrez gene",
                      "n_genes_common_all_modeling_platforms": idv["n_common_genes"]},
    "sex_check": {"method": "Y-gene score (RPS4Y1,DDX3Y,KDM5D,USP9Y,EIF1AY,UTY,ZFY,NLGN4Y) and XIST, two-means per chip type pooled across series; GSE35978 Y-only (no XIST cluster on HuGene 1.0 ST core)",
                  "per_series": sex_summary, "flagged_samples": sex_list, "donors_with_both_sexes_in_expression": donors_both,
                  "decision_written_before_modeling": "PRIMARY analysis excludes 'affected' donors: modeling sample sex conflicts with reported sex, any sample of the donor MIXED (Y and XIST both high), expression shows both sexes across the donor's samples, or confirmed identity with a different diagnosis label. Sensitivity (required control): all-eligible-donor analysis vs without affected donors vs without matched random donors (same diagnosis x cohort composition, 20 repeats). LOW_BOTH samples (weak probes, mostly U Michigan U133 Plus 2 arrays whose XIST calls male as reported) are NOT treated as affected."},
    "donor_identity": {"method": "expression residual correlation on genotype-driven identity genes (300 genes, top-10 PCs removed per unit); mutual best hits; thresholds from cross-bank negative controls",
                       "validation": idv["validation"], "calibration": idv.get("calibration_recovery"),
                       "negative_control": idv.get("negative_control"), "z_match": Z_MATCH, "z_probable": Z_PROB, "z_split_grouping": Z_SPLIT,
                       "metadata_matching": meta_m["pairs"],
                       "dedup_summary": wide,
                       "confirmed_matches_with_conflicting_diagnosis": conflicts,
                       "within_series_hits_not_explained_by_metadata": within_unexplained,
                       "confidence_statement": "Detection sensitivity is limited (validation 0.42 at probable tier for cross-tissue pairs), so the number of duplicates found is a LOWER BOUND. Splitting therefore uses liberal split groups and same-bank series are merged into one cohort."},
    "cohorts": {"Stanley": ["GSE35978 (model: parietal)", "GSE12649 (audit)"], "Pittsburgh": ["GSE53987 (model: PFC)", "GSE54567/GSE54568 (model: BA9)", "GSE54571/GSE54572 (audit: BA25)"],
                "CharingCross": ["GSE17612"], "Victoria": ["GSE21138"], "Pritzker": ["GSE92538"]},
    "region_decision": "one cortical sample per donor; prefrontal preferred; Stanley modelled from GSE35978 parietal cortex (single platform/region for the cohort). Cerebellum, hippocampus, striatum, amygdala/ACC not used for modeling. Technical replicate arrays of the same donor+region are averaged (Pritzker; repeated GSE35978 chips).",
    "class_inclusion": {"rule": "disorder enters Task A/B iff >= 3 distinct cohorts each contribute >= %d unique donors after deduplication (also required after removing affected donors)" % MIN_DONORS,
                        "min_donors_per_contributing_cohort": MIN_DONORS, "per_disorder": inclusion,
                        "LOCKED_CLASS_LIST": CLASSES, "locked": True},
    "tasks": tasks,
    "modeling_donors": {"n_eligible": int(len(model)), "n_affected": int(model.affected.sum()),
                        "affected_donors": affected_list, "unresolvable_label_donors_excluded": unresolvable,
                        "by_cohort_dx": tab(model), "by_batch": {k: int(v) for k, v in model.batch.value_counts().items()}},
}
write_json(audit, ROOT / "DATA_AUDIT.json")
print(pd.DataFrame(affected_list).to_string()); print("unresolvable:", unresolvable)
print(json.dumps({"classes": CLASSES, "tasks": {k: (v["n_total"], v["donors"]) for k, v in tasks.items()},
                  "dedup": wide, "conflicts": len(conflicts), "affected": int(model.affected.sum())}, indent=1, default=str))
