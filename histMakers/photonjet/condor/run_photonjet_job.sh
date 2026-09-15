#!/usr/bin/env bash
# One original-tree job; Condor allocates one CPU.
set -e
source "${6:-${PJ_SETUP_SCRIPT:-/sphenix/u/shuhang98/setup.sh}}"
set -uo pipefail
WORKDIR=${5:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"}
cd "$WORKDIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec python3 scripts/workflow.py worker "${1:?missing config}" --product "${2:?missing product}" --input "${3:?missing input}" --tag "${4:?missing tag}" --threads "${7:-1}"
