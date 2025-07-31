#!/bin/bash

set -euo pipefail

log() {
    echo -e "\n[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# CONFIGURATION
VFDB_DMND="vfdb/VFDB_prot.dmnd"
TADB3_DMND_DIR="tadb3/tadb3_combined.fasta"
EGGNOG_DATA_DIR="dbs/eggnog_data"
INTERPROSCAN_BIN="dbs/interproscan-5.72-103.0/interproscan.sh"
THREADS_DIAMOND=64
THREADS_EGGNOG=80

# Input FASTA directories
declare -A FASTA_DIRS=(
    ["dieguez"]="dieguez_uniref_mapped/fastas"
    ["ev"]="ev_uniref_mapped/fastas"
)

# Output directories for each tool
declare -A VFDB_OUT=(
    ["dieguez"]="vfdb/dieguez"
    ["ev"]="vfdb/ev"
)
declare -A TADB3_OUT=(
    ["dieguez"]="tadb3/dieguez"
    ["ev"]="tadb3/ev"
)
declare -A EGGNOG_OUT=(
    ["dieguez"]="eggnog/dieguez"
    ["ev"]="eggnog/ev"
)
declare -A INTERPROSCAN_OUT=(
    ["dieguez"]="interproscan/dieguez"
    ["ev"]="interproscan/ev"
)

# SETUP
log "Creating all output directories..."
for dir in "${VFDB_OUT[@]}" "${TADB3_OUT[@]}" "${EGGNOG_OUT[@]}" "${INTERPROSCAN_OUT[@]}"; do
    mkdir -p "$dir"
done

# FUNCTIONS
run_diamond_blastp() {
    local fasta_dir=$1
    local output_dir=$2
    local db=$3
    local label=$4

    for fasta in "$fasta_dir"/*.fasta; do
        local base=$(basename "$fasta" .fasta)
        local out="$output_dir/${base}_vs_$(basename "$db" .dmnd).tsv"
        log "DIAMOND VFDB [$label]: $base vs $(basename "$db")"
        diamond blastp -q "$fasta" -d "$db" -o "$out" --outfmt 6 --threads "$THREADS_DIAMOND"
    done
}

run_diamond_tadb3_all() {
    local fasta_dir=$1
    local output_dir=$2

    for fasta in "$fasta_dir"/*.fasta; do
        local base=$(basename "$fasta" .fasta)
        for dmnd in "$TADB3_DMND_DIR"/*.dmnd; do
            local dbname=$(basename "$dmnd" .dmnd)
            local out="$output_dir/${base}_vs_${dbname}.tsv"
            log "DIAMOND TADB3: $base vs $dbname"
            diamond blastp -q "$fasta" -d "$dmnd" -o "$out" --outfmt 6 --threads "$THREADS_DIAMOND"
        done
    done
}

run_eggnog_mapper() {
    local fasta_dir=$1
    local output_dir=$2

    for fasta in "$fasta_dir"/*.fasta; do
        local base=$(basename "$fasta" .fasta)
        local out="$output_dir/${base}_eggnog"
        log "EggNOG-mapper: $base"
        emapper.py -i "$fasta" -o "$out" \
            --pident 50 --query_cover 80 \
            --data_dir "$EGGNOG_DATA_DIR" \
            -m diamond --dmnd_db "$EGGNOG_DATA_DIR/eggnog_proteins.dmnd" \
            --cpu "$THREADS_EGGNOG"
    done
}

run_interproscan() {
    local fasta_dir=$1
    local output_dir=$2

    for fasta in "$fasta_dir"/*.fasta; do
        local base=$(basename "$fasta" .fasta)
        local out="$output_dir/${base}.tsv"
        log "InterProScan: $base"
        "$INTERPROSCAN_BIN" -i "$fasta" -t p -f tsv -o "$out" -pa -iprlookup -goterms
    done
}

# PIPELINE EXECUTION
log "Starting DIAMOND VFDB analysis..."
for dataset in "${!FASTA_DIRS[@]}"; do
    run_diamond_blastp "${FASTA_DIRS[$dataset]}" "${VFDB_OUT[$dataset]}" "$VFDB_DMND" "$dataset"
done
log "Completed DIAMOND VFDB analysis."

log "Starting DIAMOND TADB3 analysis..."
for dataset in "${!FASTA_DIRS[@]}"; do
    run_diamond_tadb3_all "${FASTA_DIRS[$dataset]}" "${TADB3_OUT[$dataset]}"
done
log "Completed DIAMOND TADB3 analysis."

log "Starting EggNOG-mapper analysis..."
for dataset in "${!FASTA_DIRS[@]}"; do
    run_eggnog_mapper "${FASTA_DIRS[$dataset]}" "${EGGNOG_OUT[$dataset]}"
done
log "Completed EggNOG-mapper analysis."

log "Starting InterProScan analysis..."
for dataset in "${!FASTA_DIRS[@]}"; do
    run_interproscan "${FASTA_DIRS[$dataset]}" "${INTERPROSCAN_OUT[$dataset]}"
done
log "Completed InterProScan analysis."

log "Functional annotation pipeline complete. All steps finished."

