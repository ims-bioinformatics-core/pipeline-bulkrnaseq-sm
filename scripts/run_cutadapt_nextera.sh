#!/bin/bash
set -euo pipefail

# Usage:
# run_cutadapt_nextera.sh <input_r1> <input_r2> <output_r1> <output_r2>

IN1=$1
IN2=$2
OUT1=$3
OUT2=$4

if [ -z "$IN1" ] || [ -z "$IN2" ] || [ -z "$OUT1" ] || [ -z "$OUT2" ]; then
    echo "Usage: $0 <input_r1> <input_r2> <output_r1> <output_r2>"
    exit 1
fi

mkdir -p "$(dirname "$OUT1")"
mkdir -p "$(dirname "$OUT2")"

echo "✂️ Trimming Nextera adapters"
echo "  R1: $IN1 -> $OUT1"
echo "  R2: $IN2 -> $OUT2"

cutadapt -j 0 \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -m 25 \
    --report=minimal \
    -o "$OUT1" \
    -p "$OUT2" \
    "$IN1" "$IN2"

echo "✅ Nextera adapter trimming complete"
echo "  Output R1: $OUT1"
echo "  Output R2: $OUT2"

