#!/usr/bin/env bash
# Run one histogram macro on one chunk of parts. Used by Condor and by submit_photonjet.sh --local.
#   run_photonjet_job.sh <photonjet dir> <config> <spectrum|events> <product> <chunk list> <tag>
set -e
source /opt/sphenix/core/bin/sphenix_setup.sh -n new >/dev/null 2>&1
# yaml-cpp, built outside the sPHENIX release: headers for the macros and the library for the symbol lookup.
export ROOT_INCLUDE_PATH=/sphenix/u/shuhang98/install/include:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=/sphenix/u/shuhang98/install/lib64:$LD_LIBRARY_PATH
set -uo pipefail
cd "${1:?missing photonjet directory}"
CONFIG=${2:?missing config}
CORE=${3:?missing core}
PRODUCT=${4:?missing product}
LIST=${5:?missing chunk list}
TAG=${6:?missing tag}
case "$CORE" in
    spectrum) MACRO=RecoEffCalculator_TTreeReader ;;
    events) MACRO=EvtCharacter ;;
    *) echo "unknown core $CORE" >&2; exit 2 ;;
esac
echo "[run_photonjet_job] $MACRO config=$CONFIG product=$PRODUCT list=$LIST tag=$TAG"
exec root -l -b -q "histmakers/${MACRO}.C(\"${CONFIG}\",\"${PRODUCT}\",\"${LIST}\",\"${TAG}\")"
