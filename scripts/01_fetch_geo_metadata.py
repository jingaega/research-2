"""Fetch metadata-only (no expression values) GEO records for candidate series.

Uses GEO acc.cgi with targ=all&view=brief&form=text, which returns SERIES, PLATFORM
and SAMPLE annotation blocks without data tables.
"""
import sys, time, pathlib, requests

OUT = pathlib.Path(__file__).resolve().parents[1] / "data" / "geo_meta"
OUT.mkdir(parents=True, exist_ok=True)
URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=text&view=brief"

def fetch(acc):
    f = OUT / f"{acc}.txt"
    if f.exists() and f.stat().st_size > 1000:
        return f
    for attempt in range(4):
        try:
            r = requests.get(URL.format(acc=acc), timeout=120)
            r.raise_for_status()
            f.write_text(r.text)
            return f
        except Exception as e:  # network retry
            print(acc, "retry", attempt, e, file=sys.stderr)
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError(acc)

if __name__ == "__main__":
    for acc in sys.argv[1:]:
        p = fetch(acc)
        print(acc, p.stat().st_size)
        time.sleep(0.4)
