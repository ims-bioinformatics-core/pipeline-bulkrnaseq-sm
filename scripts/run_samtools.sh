#!/bin/bash
set -euo pipefail

# ==========================
# Usage:
# ./dedup_sorted_bam.sh input.bam output_prefix threads annotation.gtf strand
#
# Example:
# ./dedup_sorted_bam.sh \
#   UDI0079rna.Aligned.sortedByCoord.out.bam \
#   UDI0079rna \
#   8 \
#   genes.gtf \
#   0
#
# Strand:
# 0 = unstranded (most common)
# 1 = stranded forward
# 2 = stranded reverse
# ==========================

INPUT_BAM="$1"
OUT_PREFIX="$2"
THREADS="${3:-4}"
GTF_FILE="$4"
STRAND="${5:-0}"

# -----------------------------
# Temporary files
# -----------------------------
TMP_NAME_SORT="${OUT_PREFIX}.name_sorted.bam"
TMP_FIXMATE="${OUT_PREFIX}.fixmate.bam"
TMP_COORD_SORT="${OUT_PREFIX}.fixmate.sorted.bam"
PRIMARY_BAM="${OUT_PREFIX}.primary.bam"

# -----------------------------
# Output BAM (deduplicated)
# -----------------------------
DEDUP_BAM="${OUT_PREFIX}.Aligned.sortedByCoord.out.bam"
COUNTS_FILE="${OUT_PREFIX}.ReadsPerGene.out.tab"

FLAGSTAT_DEDUP="${OUT_PREFIX}.dedup.flagstat.txt"
STATS_DEDUP="${OUT_PREFIX}.dedup.stats.txt"

# -----------------------------
# 1️⃣ Sort by read name
# -----------------------------
echo "▶ Sorting BAM by read name"
samtools sort -n -@ "$THREADS" -o "$TMP_NAME_SORT" "$INPUT_BAM"

# -----------------------------
# 2️⃣ Fix mate information
# -----------------------------
echo "▶ Running samtools fixmate"
samtools fixmate -m -@ "$THREADS" "$TMP_NAME_SORT" "$TMP_FIXMATE"

# -----------------------------
# 3️⃣ Sort by coordinate
# -----------------------------
echo "▶ Sorting fixmate BAM by coordinate"
samtools sort -@ "$THREADS" -o "$TMP_COORD_SORT" "$TMP_FIXMATE"

# -----------------------------
# 4️⃣ Remove duplicates (directly from coordinate-sorted BAM)
# -----------------------------
echo "▶ Removing duplicates"
samtools markdup -r -@ "$THREADS" "$TMP_COORD_SORT" "$DEDUP_BAM"

# -----------------------------
# 5️⃣ Index deduplicated BAM
# -----------------------------
echo "▶ Indexing deduplicated BAM"
samtools index "$DEDUP_BAM"

# -----------------------------
# 6️⃣ featureCounts (gene-level counts)
# -----------------------------
echo "▶ Running featureCounts"
featureCounts \
  -T "$THREADS" \
  -p \
  -s "$STRAND" \
  -t exon \
  -g gene_id \
  -a "$GTF_FILE" \
  -o "$OUT_PREFIX.ReadsPerGene.out.tab" \
  "$OUT_PREFIX.Aligned.sortedByCoord.out.bam"

# -----------------------------
# 7️⃣ QC metrics (primary alignments only)
# -----------------------------
echo "▶ Generating samtools QC metrics (primary alignments only)"
samtools view -@ "$THREADS" \
    -F 0x100 \
    -F 0x800 \
    -b "$DEDUP_BAM" > "$PRIMARY_BAM"

samtools flagstat "$PRIMARY_BAM" > "$FLAGSTAT_DEDUP"
samtools stats "$PRIMARY_BAM"    > "$STATS_DEDUP"

# -----------------------------
# 8️⃣ Cleanup temporary files
# -----------------------------
rm -f \
    "$TMP_NAME_SORT" \
    "$TMP_FIXMATE" \
    "$TMP_COORD_SORT" \
    "$PRIMARY_BAM"

echo "✅ Pipeline complete"
echo "   Final BAM   : $DEDUP_BAM"
echo "   Counts      : $COUNTS_FILE"
echo "   QC files    :"
echo "     - $FLAGSTAT_DEDUP"
echo "     - $STATS_DEDUP"
