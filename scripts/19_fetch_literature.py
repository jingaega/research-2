"""Literature retrieval via NCBI E-utilities: esearch (title/topic query) -> efetch abstract text.
Saves each record to data/lit/pmid_<id>.txt and logs query -> PMID in results/literature_queries.json."""
import json, sys, time, pathlib, requests
from common import DATA, RES, write_json
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = DATA / "lit"
queries = json.loads(sys.argv[1])
log_f = RES / "literature_queries.json"
log = json.loads(log_f.read_text()) if log_f.exists() else []
for q in queries:
    r = requests.get(E + "esearch.fcgi", params={"db": "pubmed", "term": q, "retmax": 3, "retmode": "json"}, timeout=60).json()
    ids = r["esearchresult"]["idlist"]
    log.append({"query": q, "pmids": ids})
    for pid in ids[:1]:
        t = requests.get(E + "efetch.fcgi", params={"db": "pubmed", "id": pid, "rettype": "abstract", "retmode": "text"}, timeout=60).text
        (OUT / f"pmid_{pid}.txt").write_text(t)
        print("=" * 100, "\nQUERY:", q, "\nPMID:", pid, "\n", t[:2600])
    time.sleep(0.4)
write_json(log, log_f)
