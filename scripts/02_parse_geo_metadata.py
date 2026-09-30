"""Parse GEO brief-text records into a flat per-sample table (raw characteristics kept)."""
import pathlib, json, re
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
META = ROOT / "data" / "geo_meta"

def parse(path):
    series, platforms, samples = {}, {}, []
    cur, kind = None, None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("^"):
            kind, _, acc = line[1:].partition(" = ")
            cur = {"_acc": acc.strip()}
            if kind == "SERIES": series = cur
            elif kind == "PLATFORM": platforms[cur["_acc"]] = cur
            elif kind == "SAMPLE": samples.append(cur)
            continue
        if not line.startswith("!") or cur is None: continue
        k, _, v = line[1:].partition(" = ")
        cur.setdefault(k, []).append(v)
    return series, platforms, samples

rows = []
for f in sorted(META.glob("GSE*.txt")):
    if f.stem in ("GSE35974", "GSE35977"): continue  # SubSeries of GSE35978; samples already there
    s, p, sm = parse(f)
    gse = s["_acc"]
    for x in sm:
        ch = {}
        for c in x.get("Sample_characteristics_ch1", []):
            k, sep, v = c.partition(":")
            key = k.strip().lower() if sep else "unlabelled"
            ch[key] = (ch[key] + " | " + v.strip()) if key in ch else v.strip()
        rows.append({"gse": gse, "gsm": x["_acc"],
                     "series_in_record": ";".join(x.get("Sample_series_id", [])),
                     "platform": x.get("Sample_platform_id", [""])[0],
                     "title": x.get("Sample_title", [""])[0],
                     "source": x.get("Sample_source_name_ch1", [""])[0],
                     "organism": x.get("Sample_organism_ch1", [""])[0],
                     "supp": ";".join(x.get("Sample_supplementary_file", [])),
                     "description": " | ".join(x.get("Sample_description", [])),
                     "char": json.dumps(ch)})
df = pd.DataFrame(rows)
df.to_csv(ROOT / "data" / "samples_raw.tsv", sep="\t", index=False)
print(df.groupby(["gse", "platform"]).size())
