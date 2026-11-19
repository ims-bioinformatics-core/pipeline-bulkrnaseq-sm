#!/usr/bin/env bash

# Usage:
# ./utils-downsample_fastq.sh <fastq_dir> <target_reads> [output_dir]
# Example:
# ./utils-downsample_fastq.sh ./raw_data 20000000 downsampled_reads

FASTQ_DIR="$1"
TARGET_READS="$2"
OUT_DIR="${3:-downsampled}"  # default output directory

# ---------------------------
# Input checks
# ---------------------------
if [[ -z "$FASTQ_DIR" || -z "$TARGET_READS" ]]; then
    echo "Usage: $0 <fastq_dir> <target_reads> [output_dir]"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: $FASTQ_DIR is not a directory"
    exit 1
fi

# Create output directory
mkdir -p "$OUT_DIR"

echo "Processing FASTQ files in $FASTQ_DIR ..."
echo "Target reads per sample: $TARGET_READS"
echo "Output directory: $OUT_DIR"

# ---------------------------
# Process each sample
# ---------------------------
for fastq1 in "$FASTQ_DIR"/*.r_1.fq.gz; do
    sample=$(basename "$fastq1" .r_1.fq.gz)
    fastq2="$FASTQ_DIR/${sample}.r_2.fq.gz"

    if [[ ! -f "$fastq2" ]]; then
        echo "Warning: Paired FASTQ file for sample $sample not found, skipping."
        continue
    fi

    # Count total reads using seqtk (fast)
    total_reads=$(seqtk seq -l0 "$fastq1" | grep -c "^@")

    if (( total_reads <= TARGET_READS )); then
        echo "Sample $sample: total reads ($total_reads) <= target ($TARGET_READS), skipping downsample."
        continue
    fi

    # Calculate fraction for seqtk
    frac=$(awk -v t="$TARGET_READS" -v tot="$total_reads" 'BEGIN{printf "%.6f", t/tot}')

    echo "Downsampling $sample: total reads $total_reads -> target $TARGET_READS (~$frac fraction)"

    out1="$OUT_DIR/${sample}.r_1.fq.gz"
    out2="$OUT_DIR/${sample}.r_2.fq.gz"

    seqtk sample -s100 "$fastq1" "$frac" | gzip > "$out1"
    seqtk sample -s100 "$fastq2" "$frac" | gzip > "$out2"

    echo "Downsampled FASTQ files for $sample created: $out1, $out2"
done

echo "✅ Done."
