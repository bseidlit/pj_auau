#!/usr/bin/env bash
# Write lists/<system>_<product>.list from the dated v21 exports, one path per line.
# A part is identified by the number in its file name, not by its place in the list.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE=/sphenix/tg/tg01/bulk/jbennett/PhotonJetTrees_v2

# Products the producer has promoted to a top-level parts/ directory ship a files.txt, which is
# the complete plain path list. The others are still dated exports, listed here.
declare -A EXPORTS=(
    [pp_data]=data/pp
    [pp_sim_signal]=simulation/pp_photonjet
    [pp_sim_inclusive]=simulation/pp_inclusive
    [pp_sim_di]=simulation/pp_di
    [auau_data]=data/auau
    [auau_sim_signal]=simulation/auau_embedded_photonjet
    [auau_sim_inclusive]=simulation/auau_embedded_inclusive
)

mkdir -p "$HERE/lists"
for name in $(printf '%s\n' "${!EXPORTS[@]}" | sort); do
    product="$RELEASE/${EXPORTS[$name]}"
    if [[ -f $product/files.txt ]]; then
        cp "$product/files.txt" "$HERE/lists/$name.list"
    elif compgen -G "$product/parts/part_*.root" > /dev/null; then
        find "$product/parts" -maxdepth 1 -name 'part_*.root' | sort > "$HERE/lists/$name.list"
    else
        echo "skip $name: nothing in $product" >&2
        continue
    fi
    echo "$name: $(wc -l < "$HERE/lists/$name.list") parts"
done
# The pp double-interaction parts mix photon and jet samples, so they join both simulation lists
# and each macro keeps the events of its product's kind of sample (Job::OwnsSample).
if [[ -s $HERE/lists/pp_sim_di.list ]]; then
    for name in pp_sim_signal pp_sim_inclusive; do
        cat "$HERE/lists/pp_sim_di.list" >> "$HERE/lists/$name.list"
        echo "$name: $(wc -l < "$HERE/lists/$name.list") parts with the DI parts"
    done
fi
