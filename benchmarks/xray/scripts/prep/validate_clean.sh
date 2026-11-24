#!/bin/bash
# Run the pipeline's own cleaning in strict mode (as run-target --strict-clean does on the cluster) on every target.
# Usage: validate_clean.sh <inputs dir> <scratch out dir> [pipeline CLI]
set -uo pipefail
inputs=$1; out=$2; cli=${3:-af2rank-pipeline}
mkdir -p "$out"; ok=0; failed=()
for dir in "$inputs"/*/; do
  t=$(basename "$dir"); rm -rf "$out/$t"
  if "$cli" clean --target "$t" --fasta "$dir/target.fasta" --models "$dir/models" --out "$out/$t" --strict-clean > "$out/$t.log" 2>&1; then
    ok=$((ok + 1))
  else
    failed+=("$t"); echo "FAILED $t: $(tail -1 "$out/$t.log")"
  fi
done
echo "strict clean: $ok passed, ${#failed[@]} failed: ${failed[*]:-none}"
