#!/usr/bin/env bash
# Direct original-tree workflow. Extra arguments are forwarded to each step.
# Source the analysis environment first. Each worker uses one CPU.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG=${1:?usage: oneforall_photonjet.sh <config> <chunks|local|submit|status|merge|check|yield> [options]}
STEP=${2:?missing step}
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
shift 2
cd "$HERE"
case "$STEP" in
  prepare) echo "prepare is retired; run chunks to read original ROOT inputs directly" >&2; exit 2 ;;
  chunks) exec python3 "$HERE/scripts/make_chunks.py" "$CONFIG" "$@" ;;
  local|merge|check|yield|status) exec python3 "$HERE/scripts/workflow.py" "$STEP" "$CONFIG" "$@" ;;
  submit)
    if [[ ${PJ_THREADS:-1} != 1 ]]; then
      echo "the direct maker requires one CPU; use Condor jobs for parallelism" >&2; exit 2
    fi
    VAR=$(python3 "$HERE/scripts/cfgval.py" "$CONFIG" output.var_type)
    CHUNKS=${PJ_CHUNK_DIR:-$HERE/chunks}
    LOGS=${PJ_LOG_DIR:-$HERE/logs}
    test -f "$CHUNKS/$VAR/jobs.json"
    mkdir -p "$LOGS"
    cd "$HERE"
    exec condor_submit -append "photonjet_dir=$HERE" -append "setup_script=${PJ_SETUP_SCRIPT:-/sphenix/u/shuhang98/setup.sh}" -append "list_file=$CHUNKS/$VAR/jobs_all.list" -append "log_dir=$LOGS" condor/submit_photonjet.sub "$@" ;;
  *) echo "unknown step: $STEP" >&2; exit 1 ;;
esac
