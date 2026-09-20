#!/usr/bin/env bash
# Split one product's part list into chunks and run one job per chunk, on Condor or here.
#
#   ./submit_photonjet.sh <config.yaml> <spectrum|events> <data|sim_signal|sim_inclusive> [parts per job] [--local]
#
# Each job writes <output.directory>/<core>_<product>_<var_type>.chunkNNNNN.root.
# Merge them with merge_photonjet.sh once every job has finished.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/photonjet_common.sh"

[[ $# -ge 3 ]] || { sed -n '2,8p' "$0" >&2; exit 2; }
read_job "$1" "$2" "$3"
PARTS_PER_JOB=50
LOCAL=0
for option in "${@:4}"; do
    if [[ $option == --local ]]; then LOCAL=1; else PARTS_PER_JOB=$option; fi
done
[[ $PARTS_PER_JOB =~ ^[1-9][0-9]*$ ]] || { echo "parts per job must be a positive integer" >&2; exit 2; }

[[ -s $MASTER_LIST ]] || { echo "missing part list $MASTER_LIST (run make_lists.sh)" >&2; exit 1; }
[[ ! -e $CHUNK_DIR ]] || { echo "$CHUNK_DIR exists: this product was already submitted; remove it to redo" >&2; exit 1; }
if compgen -G "$OUTPUT_STEM.*" > /dev/null; then
    echo "outputs already exist for $OUTPUT_STEM; use a new var_type" >&2; exit 1
fi

mkdir -p "$CHUNK_DIR" "$HERE/logs"
split -l "$PARTS_PER_JOB" -d -a 5 --additional-suffix=.list "$MASTER_LIST" "$CHUNK_DIR/chunk"
(cd "$CHUNK_DIR" && ls chunk*.list | sed 's/\.list$//') > "$CHUNK_DIR/jobs.txt"
echo "$(wc -l < "$CHUNK_DIR/jobs.txt") jobs of up to $PARTS_PER_JOB parts: $CHUNK_DIR"

if [[ $LOCAL == 1 ]]; then
    while read -r chunk; do
        "$HERE/condor/run_photonjet_job.sh" "$HERE" "$CONFIG" "$CORE" "$PRODUCT" "$CHUNK_DIR/$chunk.list" "$chunk"
    done < "$CHUNK_DIR/jobs.txt"
else
    condor_submit -append "photonjet_dir=$HERE" -append "config=$CONFIG" -append "core=$CORE" \
        -append "product=$PRODUCT" -append "chunk_dir=$CHUNK_DIR" "$HERE/condor/submit_photonjet.sub"
fi
