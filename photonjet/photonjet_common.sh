# Shared by submit_photonjet.sh and merge_photonjet.sh: resolve one (config, core, product) job.
# Sets CONFIG CORE PRODUCT MASTER_LIST CHUNK_DIR OUTPUT_STEM; the final file is $OUTPUT_STEM.root.

config_value() {  # config_value <config.yaml> <section> <key>
    python3 -c 'import sys
from ruamel.yaml import YAML
print(YAML(typ="safe").load(open(sys.argv[1]))[sys.argv[2]][sys.argv[3]])' "$@"
}

read_job() {  # read_job <config.yaml> <core> <product>
    CONFIG="$(readlink -e "$1")" || { echo "missing config $1" >&2; exit 1; }
    CORE=$2
    PRODUCT=$3
    [[ $CORE == spectrum || $CORE == events ]] || { echo "core must be spectrum or events" >&2; exit 2; }
    [[ $PRODUCT == data || $PRODUCT == sim_signal || $PRODUCT == sim_inclusive ]] ||
        { echo "product must be data, sim_signal or sim_inclusive" >&2; exit 2; }
    local variation directory
    variation="$(config_value "$CONFIG" output var_type)"
    directory="$(config_value "$CONFIG" output directory)"
    MASTER_LIST="$(config_value "$CONFIG" photonjet "${PRODUCT}_files_list")"
    # Relative paths in the config are relative to this directory, where the macros run.
    [[ $directory == /* ]] || directory="$HERE/$directory"
    [[ $MASTER_LIST == /* ]] || MASTER_LIST="$HERE/$MASTER_LIST"
    CHUNK_DIR="${PHOTONJET_CHUNKS:-$HERE/chunks}/$variation/${CORE}_$PRODUCT"
    OUTPUT_STEM="$directory/${CORE}_${PRODUCT}_$variation"
}
