#!/bin/bash
# merge_lanes.sh
# Usage: bash merge_lanes.sh /path/to/fastq_dir

set -euo pipefail

INPUT_DIR="$1"
OUTPUT_DIR="${INPUT_DIR%/}/merged_lanes"

mkdir -p "$OUTPUT_DIR"

echo "🔍 Searching for FASTQ files in: $INPUT_DIR"
echo "📁 Merged files will be saved in: $OUTPUT_DIR"

# Extract unique sample names (everything before third underscore)
samples=$(find "$INPUT_DIR" -type f -name "*.fq.gz" \
    | xargs -n1 basename \
    | awk -F'_' '{print $1"_"$2}' \
    | sort -u)

for sample in $samples; do
    echo "🔄 Processing sample: $sample"

    # Merge R1 files (strip lane & barcode part)
    R1_files=$(find "$INPUT_DIR" -type f -name "${sample}_*_*_1.fq.gz" | sort)
    if [ -n "$R1_files" ]; then
        echo "   📥 Merging R1 files..."
        cat $R1_files > "$OUTPUT_DIR/${sample}_1.fq.gz"
    fi

    # Merge R2 files (strip lane & barcode part)
    R2_files=$(find "$INPUT_DIR" -type f -name "${sample}_*_*_2.fq.gz" | sort)
    if [ -n "$R2_files" ]; then
        echo "   📥 Merging R2 files..."
        cat $R2_files > "$OUTPUT_DIR/${sample}_2.fq.gz"
    fi

    echo "✅ Merged files saved for sample: $sample"
done

echo "🎯 All merging complete."
