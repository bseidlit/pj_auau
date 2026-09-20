#!/usr/bin/env bash
# Merge one product's chunk outputs into <output.directory>/<core>_<product>_<var_type>.root.
#
#   ./merge_photonjet.sh <config.yaml> <spectrum|events> <data|sim_signal|sim_inclusive>
#
# Refuses unless every submitted chunk wrote its file, so a failed or still-running job
# cannot produce a silently partial merge. FinalizeMerged.C then checks the event count
# and builds the efficiencies and the response once, from the complete counts.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/photonjet_common.sh"

[[ $# -eq 3 ]] || { sed -n '2,8p' "$0" >&2; exit 2; }
read_job "$1" "$2" "$3"
[[ -s $CHUNK_DIR/jobs.txt ]] || { echo "no submitted chunks in $CHUNK_DIR" >&2; exit 1; }
[[ ! -e $OUTPUT_STEM.root ]] || { echo "merged output already exists: $OUTPUT_STEM.root" >&2; exit 1; }

missing=0
: > "$CHUNK_DIR/outputs.txt"
while read -r chunk; do
    output="$OUTPUT_STEM.$chunk.root"
    if [[ -s $output ]]; then echo "$output" >> "$CHUNK_DIR/outputs.txt"; else echo "missing $output" >&2; missing=$((missing + 1)); fi
done < "$CHUNK_DIR/jobs.txt"
[[ $missing -eq 0 ]] || { echo "$missing of $(wc -l < "$CHUNK_DIR/jobs.txt") chunks are missing; not merging" >&2; exit 1; }

temporary="$OUTPUT_STEM.merging.$$.root"
trap 'rm -f "$temporary"' EXIT
hadd "$temporary" @"$CHUNK_DIR/outputs.txt" > "$CHUNK_DIR/hadd.log" 2>&1 || { tail -5 "$CHUNK_DIR/hadd.log" >&2; exit 1; }
(cd "$HERE" && root -l -b -q "histmakers/FinalizeMerged.C(\"$temporary\",\"$CHUNK_DIR/outputs.txt\")")
mv -n "$temporary" "$OUTPUT_STEM.root"
[[ ! -e $temporary ]] || { echo "merged output appeared meanwhile: $OUTPUT_STEM.root" >&2; exit 1; }
echo "merged $(wc -l < "$CHUNK_DIR/outputs.txt") chunks: $OUTPUT_STEM.root"
