#!/bin/bash

# Check for input arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <output_csv_filename>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
CSV_FILENAME="$3"

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Directory where the script is executed
RUN_DIR=$(pwd)

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Output CSV file path
OUTPUT_FILE="$OUTPUT_DIR/$CSV_FILENAME"
TMP_FILE="${OUTPUT_FILE}.tmp"

# Create temporary file
> "$TMP_FILE"

echo "Searching FASTQ files in: $INPUT_DIR"
echo "Saving relative paths from: $RUN_DIR"

# Find all R1 files from both naming conventions
find "$INPUT_DIR" -type f \( -name "*.r_1.fq.gz" -o -name "*_1.fq.gz" \) | while read -r R1; do

    BASENAME=$(basename "$R1")

    if [[ "$BASENAME" == *.r_1.fq.gz ]]; then
        # Old format:
        # example.sample.r_1.fq.gz
        R2="${R1/.r_1.fq.gz/.r_2.fq.gz}"

        SAMPLE=$(echo "$BASENAME" | cut -d'.' -f2)

    elif [[ "$BASENAME" =~ _1\.fq\.gz$ ]]; then
        # New format:
        # barcode_sample_1.fq.gz
        R2="${R1%_1.fq.gz}_2.fq.gz"

        # Barcode = first two underscore-separated fields
        SAMPLE=$(echo "$BASENAME" | awk -F'_' '{print $1"_"$2}')

    else
        continue
    fi

    # Check matching R2 exists
    if [ -f "$R2" ]; then

        # Save paths relative to current working directory
        R1_PATH=$(realpath --relative-to="$RUN_DIR" "$R1")
        R2_PATH=$(realpath --relative-to="$RUN_DIR" "$R2")

        echo "$SAMPLE,$R1_PATH,$R2_PATH" >> "$TMP_FILE"

    else
        echo "Warning: Missing R2 file for $R1" >&2
    fi

done

# Write header and sorted contents
{
    echo "Barcode,fastq_1,fastq_2"
    sort -t',' -k1,1 -k2,2 "$TMP_FILE"
} > "$OUTPUT_FILE"

rm -f "$TMP_FILE"

echo "CSV file generated: $OUTPUT_FILE"