"""Standardise per-sample metadata (diagnosis, region, sex, age, PMI, pH, within-series donor key)
from GEO characteristics only. Output: data/samples.tsv"""
import json, re
import numpy as np, pandas as pd
from common import DATA

raw = pd.read_csv(DATA / "samples_raw.tsv", sep="\t", dtype=str).fillna("")
AUDITED = ["GSE12649", "GSE17612", "GSE21138", "GSE28475", "GSE28521", "GSE29555", "GSE35978", "GSE53987",
           "GSE54567", "GSE54568", "GSE54571", "GSE54572", "GSE11223", "GSE92538", "GSE21935", "GSE5388"]
DX = {"schizophrenia": "SCZ", "scz": "SCZ", "schizophrenic": "SCZ", "schiz": "SCZ",
      "bipolar": "BD", "bipolar disorder": "BD", "major depressive disorder": "MDD", "depression": "MDD",
      "mdd case": "MDD", "mdd": "MDD", "control": "CTL", "controls": "CTL", "unaffected": "CTL",
      "healthy control": "CTL", "autism": "ASD", "alcoholic": "AUD"}

def num(x):
    m = re.search(r"-?\d+(\.\d+)?", str(x)); return float(m.group()) if m else np.nan

def sexnorm(x):
    x = str(x).strip().lower()
    return {"m": "M", "male": "M", "f": "F", "female": "F"}.get(x, "")

rows = []
for _, r in raw[raw.gse.isin(AUDITED)].iterrows():
    c = json.loads(r.char); g = r.gse
    d = dict(gse=g, gsm=r.gsm, platform=r.platform, title=r.title, dx="", region="", sex_reported="",
             age=np.nan, pmi=np.nan, ph=np.nan, donor_local="", sample_kind="tissue", note="")
    if g == "GSE11223":
        d.update(dx="NON_PSYCH", region="sigmoid colon (non-brain)", sample_kind="non_brain")
    elif g == "GSE12649":
        d["dx"] = DX[r.source.strip().lower()]; d["region"] = "PFC_BA46"
        d["donor_local"] = r.title.split("BA46-")[1]
    elif g in ("GSE17612", "GSE21935"):
        d["dx"] = "SCZ" if "Scz" in r.title or "schizo" in r.source.lower() else "CTL"
        d["region"] = "PFC_BA10" if g == "GSE17612" else "STG_BA22"
        p = r.title.split("_"); d["donor_local"] = p[0]
        d["sex_reported"] = sexnorm(c.get("gender", "")) or sexnorm(p[2])
        if not c.get("gender"): d["note"] = "sex missing in characteristics; taken from title"
        d["age"] = num(c.get("age")); d["pmi"] = num(c.get("post-mortem delay")); d["ph"] = num(c.get("ph"))
    elif g == "GSE21138":
        d["dx"] = "SCZ" if "schizophrenia" in c["stage of illness [short doi=<5 yrs; intermediate doi=7-18yrs; long doi=>28 yrs]"] else "CTL"
        d["region"] = "PFC_BA46"; d["donor_local"] = r.title
        d.update(sex_reported=sexnorm(c["sex"]), age=num(c["age"]), pmi=num(c["pmi (hrs)"]), ph=num(c["tissue ph"]))
        if r.title == "Control-7": d["note"] = "submitter flags Cont-7 as outlier removed from publication"
    elif g in ("GSE28475", "GSE28521"):
        ds = c.get("diagnosis", c.get("disease status", "")).strip().lower()
        d["dx"] = DX.get(ds, "")
        if g == "GSE28521":
            d["region"] = {"Frontal cortex": "FC", "Temporal cortex": "TC", "Cerebellum": "CB"}[c["tissue (brain region)"]]
            d["donor_local"] = r.title.split("_")[1]
        else:
            d["region"] = "brain_unspecified"
            if not ds: d.update(sample_kind="reference_rna")
            elif "Fixed" in r.title: d["sample_kind"] = "fixed_tissue"
            m = re.search(r"_(\d+[A-C]|EC\d|M\d+)", r.title); d["donor_local"] = m.group(1) if m else ""
    elif g == "GSE29555":
        d["dx"] = "AUD" if "Alcoholic" in r.source else "CTL"; d["region"] = c["tissue"]
        d["donor_local"] = re.match(r"([AC]\d+)R", r.title).group(1)
        d.update(sex_reported=sexnorm(c["gender"]), age=num(c["age"]))
    elif g == "GSE35978":
        if r.title.startswith("UNIVREP"): d.update(sample_kind="reference_rna", dx="REF")
        else:
            ds = c["disease status"].strip().lower()
            d["dx"] = {"bipolar (not bipolar)": "AMBIG", "na": "NA"}.get(ds, DX.get(ds, ds))
        d["region"] = {"parietal cortex": "PARIETAL", "cerebellum": "CB"}[c["tissue"]]
        d.update(sex_reported=sexnorm(c["sex"]), age=num(c["age"]), pmi=num(c["post-mortem interval (pmi)"]), ph=num(c["ph"]))
    elif g == "GSE53987":
        d["dx"] = DX[c["disease state"]]
        d["region"] = {"Pre-frontal cortex (BA46)": "PFC_BA46", "hippocampus": "HPC", "Associative striatum": "STR"}[c["tissue"]]
        d.update(sex_reported=sexnorm(c["gender"]), age=num(c["age"]), pmi=num(c["pmi"]), ph=num(c["ph"]))
        d["race"] = c["race"]; d["rin"] = num(c["rin"])
    elif g in ("GSE54567", "GSE54568", "GSE54571", "GSE54572"):
        m = re.search(r"(BA9|BA25)_([MF]) #(\d+)", r.title)
        d["dx"] = "MDD" if "MDD" in c["disease state"] else "CTL"
        d["region"] = "PFC_BA9" if m.group(1) == "BA9" else "ACC_BA25"
        d["sex_reported"] = m.group(2)  # sex only encoded in the series split / title
        d["donor_local"] = f"{m.group(2)}pair{m.group(3)}_{d['dx']}"
        d["note"] = "sex from series title (BA*_M/BA*_F); description text mentions temporal pole (apparent copy error)"
    elif g == "GSE92538":
        d["dx"] = DX[c["diagnosis"].lower()] if c["diagnosis"] != "Bipolar Disorder" else "BD"
        d["region"] = "PFC_DLPFC"; d["donor_local"] = c["subject id"]
        d.update(sex_reported=sexnorm(c["gender"]), age=num(c["age"]), pmi=num(c["post-mortem interval"]),
                 ph=num(c.get("tissue ph (cerebellum)")))
        d["chip"] = "U133Plus2" if r.platform == "GPL10526" else "U133A"
        d["site"] = c["site of processing"]
    elif g == "GSE5388":
        d["dx"] = "BD" if "Bipolar" in c["disease_status"] else "CTL"; d["region"] = "PFC_DLPFC"
        d.update(sex_reported=sexnorm(c["gender"]), age=num(c["age (years)"]),
                 pmi=num(c["post mortem interval (hours)"]), ph=num(c["brain ph"]))
        d["donor_local"] = r.title.split("_")[-1]
    rows.append(d)

s = pd.DataFrame(rows)
# GSE35978 / GSE53987: no donor id in GEO -> donor key from identical (dx, sex, age, pmi, pH[, race]) tuples
for g, cols in {"GSE35978": ["dx", "sex_reported", "age", "pmi", "ph"],
                "GSE53987": ["dx", "sex_reported", "age", "pmi", "ph", "race"]}.items():
    m = (s.gse == g) & (s.sample_kind == "tissue") & ~s.dx.isin(["NA"])
    key = s.loc[m, cols].astype(str).agg("|".join, axis=1)
    s.loc[m, "donor_local"] = "k" + pd.factorize(key)[0].astype(str)
na = (s.gse == "GSE35978") & (s.dx == "NA")  # no metadata at all: one provisional donor per sample
s.loc[na, "donor_local"] = "na_" + s.loc[na, "gsm"]
s.to_csv(DATA / "samples.tsv", sep="\t", index=False)
t = s[s.sample_kind == "tissue"]
print(t.groupby(["gse", "region", "dx"]).agg(samples=("gsm", "size"), donors=("donor_local", "nunique")).to_string())
for g in ["GSE35978", "GSE53987"]:
    x = t[t.gse == g].groupby("donor_local").agg(n=("gsm", "size"), regions=("region", lambda v: ",".join(sorted(v))))
    print(g, "samples per metadata-donor:", x.n.value_counts().to_dict(), "| dup region within donor:",
          int((x.regions.str.split(",").apply(lambda v: len(v) != len(set(v)))).sum()))
