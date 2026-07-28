#!/bin/bash
set -euo pipefail

# Usage:
# ./bam_stats.sh input_dedup.bam threads
# Example:
# ./bam_stats.sh sample.dedup.bam 4

DEDUP_BAM="$1"
THREADS="${2:-4}"

# Output files in the same directory as BAM
OUT_PREFIX="${DEDUP_BAM%.bam}"  # strip .bam
PRIMARY_BAM="${OUT_PREFIX}.primary.bam"
FLAGSTAT="${OUT_PREFIX}.dedup.flagstat.txt"
STATS="${OUT_PREFIX}.dedup.stats.txt"

echo "▶ Generating stats for sample: $(basename "$OUT_PREFIX")"
echo "  Input BAM : $DEDUP_BAM"
echo "  Output prefix : $OUT_PREFIX"

echo "▶ Extracting primary alignments from BAM"
samtools view -@ "$THREADS" -F 0x100 -F 0x800 -b "$DEDUP_BAM" > "$PRIMARY_BAM"

echo "▶ Running samtools flagstat"
samtools flagstat "$PRIMARY_BAM" > "$FLAGSTAT"

echo "▶ Running samtools stats"
samtools stats "$PRIMARY_BAM" > "$STATS"

# Cleanup
rm -f "$PRIMARY_BAM"

echo "✅ Stats complete"
echo "  Flagstat : $FLAGSTAT"
echo "  Stats    : $STATS"
