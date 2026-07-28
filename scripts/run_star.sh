#!/bin/bash
set -euo pipefail


# --- Inputs ---
THREADS="$1"
GENOME_DIR=$(realpath "$2")
R1=$(realpath "$3")
R2=$(realpath "$4")
OUT_DIR=$(realpath "$5")
OUT_PREFIX="$6"
SAM_FORMAT="$7"
SAM_SORT="$8"
STRANDFIELD="$9"
FILTER_SCORE="${10}"
FILTER_MATCH="${11}"
QUANT_MODE="${12}"

# --- Prepare output directories ---
#mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

ulimit -n 20000

echo "🚀 STAR run:"
echo "    Threads: $THREADS"
echo "    Genome: $GENOME_DIR"
echo "    R1: $R1"
echo "    R2: $R2"
echo "    Output dir: $OUT_DIR"
echo "    Output prefix: $OUT_PREFIX"
echo "    SAM format: $SAM_FORMAT"
echo "    SAM sort: $SAM_SORT"
echo "    Strand field: $STRANDFIELD"
echo "    Filter Score Min Over Lread: $FILTER_SCORE"
echo "    Filter Match N Min Over Lread: $FILTER_MATCH"
echo "    Quant mode: $QUANT_MODE"

# --- Run STAR ---
STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$OUT_PREFIX." \
    --outTmpDir "$OUT_DIR/$OUT_PREFIX" \
    --outSAMtype "$SAM_FORMAT" "$SAM_SORT" \
    --outSAMstrandField "$STRANDFIELD" \
    --outFilterScoreMinOverLread "$FILTER_SCORE" \
    --outFilterMatchNminOverLread "$FILTER_MATCH" \
    --quantMode "$QUANT_MODE"
