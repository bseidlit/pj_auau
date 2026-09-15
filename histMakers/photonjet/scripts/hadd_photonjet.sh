#!/usr/bin/env bash
# Checked original-input merge, MBD import, and generation publication.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG=${1:?usage: hadd_photonjet.sh <config> [products]}
exec python3 "$HERE/workflow.py" merge "$CONFIG" --products "${2:-data,sim_signal,sim_inclusive}"
