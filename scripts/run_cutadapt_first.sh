#!/bin/bash
set -euo pipefail

# Usage:
# trim_fixed_bp.sh <input_r1> <input_r2> <trim_bp> <mode> <output_r1> <output_r2>
#
# mode = r1 | r2 | both

IN1=$1
IN2=$2
TRIM_BP=${3:-3}
MODE=${4:-r1}
OUT1=$5
OUT2=$6

if [ -z "$IN1" ] || [ -z "$OUT1" ]; then
    echo "Usage: $0 <input_r1> <input_r2> <trim_bp> <mode> <output_r1> <output_r2>"
    exit 1
fi

mkdir -p "$(dirname "$OUT1")"
mkdir -p "$(dirname "$OUT2")"

echo "✂️ Fixed-base trimming"
echo "    Trim BP: $TRIM_BP"
echo "    Mode: $MODE"
echo "    R1: $IN1"
echo "    R2: ${IN2:-NONE}"

case "$MODE" in

    r1)
        echo "Trimming R1 only"
        cutadapt -j 0 -u "$TRIM_BP" -m 25 --report=minimal \
            -o "$OUT1" "$IN1"

        [[ -n "${IN2:-}" && -f "$IN2" ]] && cp "$IN2" "$OUT2"
        ;;

    r2)
        if [[ -z "${IN2:-}" || ! -f "$IN2" ]]; then
            echo "❌ R2 requested but not provided"
            exit 1
        fi

        echo "Trimming R2 only"
        cutadapt -j 0 -u "$TRIM_BP" -m 25 --report=minimal \
            -o "$OUT2" "$IN2"

        cp "$IN1" "$OUT1"
        ;;

    both)
        if [[ -z "${IN2:-}" || ! -f "$IN2" ]]; then
            echo "❌ Both requested but R2 missing"
            exit 1
        fi

        echo "Trimming both reads"
        cutadapt -j 0 \
            -u "$TRIM_BP" \
            -U "$TRIM_BP" \
            -m 25 \
            --report=minimal \
            -o "$OUT1" \
            -p "$OUT2" \
            "$IN1" "$IN2"
        ;;

    *)
        echo "❌ Invalid mode: $MODE"
        echo "Valid modes: r1 | r2 | both"
        exit 1
        ;;

esac

echo "✅ Trimming complete"
echo "    Output R1: $OUT1"
echo "    Output R2: $OUT2"
