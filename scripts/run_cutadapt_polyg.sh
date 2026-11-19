#!/bin/bash

# Usage: trim_cutadapt.sh <input_r1> <input_r2> <output_r1> <output_r2>
# Trims poly-G tails (>=10 Gs) from paired-end reads

IN1=$1
IN2=$2
OUT1=$3
OUT2=$4

if [ -z "$IN1" ] || [ -z "$IN2" ] || [ -z "$OUT1" ] || [ -z "$OUT2" ]; then
    echo "Usage: $0 <input_r1> <input_r2> <output_r1> <output_r2>"
    exit 1
fi

# Create parent directories
mkdir -p "$(dirname "$OUT1")"
mkdir -p "$(dirname "$OUT2")"

# Poly-G trimming (3' end), removes stretches of 10 or more consecutive Gs
POLYG_ARGS=("-a" "G{10}" "-A" "G{10}")

echo "Trimming poly-G tails from paired files:"
echo "  R1: $IN1 -> $OUT1"
echo "  R2: $IN2 -> $OUT2"

# Run cutadapt
cutadapt -j 0 \
    -m 25 --length 100 --report=minimal \
    "${POLYG_ARGS[@]}" \
    -o "$OUT1" -p "$OUT2" "$IN1" "$IN2"

echo "✅ Poly-G trimming complete. Output saved as:"
echo "  R1: $OUT1"
echo "  R2: $OUT2"
