#!/bin/bash
# Sequentially download the biological assembly mmCIF for every dimer in xray_dimers.csv into ../raw/.
# Uses /usr/bin/curl (the Anaconda CA bundle fails on the RCSB certificate chain). Skips files already present.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
csv=$here/../xray_dimers.csv
raw=$here/../raw
mkdir -p "$raw"
/usr/bin/python3 -I -c '
import csv,sys
for r in csv.DictReader(open(sys.argv[1])): print(r["pdb_id"], r["assembly_id"] or "1")
' "$csv" | while read -r id asm; do
  out=$raw/${id}-assembly${asm}.cif
  [ -s "$out" ] && continue
  /usr/bin/curl -sS --fail --retry 3 -o "$out.part" "https://files.rcsb.org/download/${id}-assembly${asm}.cif" && mv "$out.part" "$out" || echo "FAILED $id"
  sleep 0.3
done
ls "$raw"/*.cif | wc -l
