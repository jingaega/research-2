#!/bin/bash
# Download GEO RAW (CEL) tarballs for series retained for modeling/QC.
cd "$(dirname "$0")/../data/raw"
for acc in "$@"; do
  stub="${acc:0:${#acc}-3}nnn"; f="${acc}_RAW.tar"
  url="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${acc}/suppl/${f}"
  if [ ! -s "$f.done" ]; then
    for i in 1 2 3 4 5; do curl -sS -C - -o "$f" "$url" && touch "$f.done" && break; sleep $((2**i)); done
  fi
  echo "$acc $(stat -c %s $f 2>/dev/null) $(date +%T)"
done
