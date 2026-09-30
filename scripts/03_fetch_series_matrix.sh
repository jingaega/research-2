#!/bin/bash
# Download deposited series-matrix files (used ONLY for label-free QC: sex check and donor identity).
cd "$(dirname "$0")/../data/matrix"
for acc in "$@"; do
  stub="${acc:0:${#acc}-3}nnn"
  url="https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${acc}/matrix/"
  for f in $(curl -sS "$url" | grep -o "${acc}[^\"]*series_matrix.txt.gz" | sort -u); do
    [ -s "$f" ] || for i in 1 2 3 4; do curl -sS -o "$f" "$url$f" && break; sleep $((2**i)); done
    echo "$f $(stat -c %s $f)"
  done
done
