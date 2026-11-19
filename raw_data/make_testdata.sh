#!/bin/bash
# Build the public mouse test data set shipped in raw_data/.
#
# Source: DDBJ BioProject PRJDB17560 (public, via ENA)
#   "Sex differences in the function of hepatic VLDLR during protein deficiency in mice"
#   Mouse liver, C57BL/6J WT vs Vldlr KO, control (20% casein) vs low-protein (2% casein) diet,
#   female and male, Illumina NovaSeq 6000, paired-end.
#
# Each run is downsampled to 500,000 read pairs with seqtk (same seed for R1 and R2 keeps
# pairs matched) and renamed to the SLX-11111.UDIxxx.AAA111.r_N.fq.gz scheme used in inputs/.
#
# Usage (from the pipeline root):  bash raw_data/make_testdata.sh
set -euo pipefail
source /rds/project/rds-O11U8YqSuCk/bioinformatics/envs/miniconda3/etc/profile.d/conda.sh
conda activate pipeline-bulkrnaseq-sm_env

OUT=raw_data
TMP=raw_data_download          # full-size FASTQs, deleted at the end (gitignored)
N=500000
SEED=11
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv&accession=PRJDB17560"
mkdir -p "$TMP" "$OUT"
curl -sf "$ENA" > "$TMP/filereport.tsv"

# barcode  run
MAP="UDI001 DRR530875
UDI002 DRR530877
UDI003 DRR530879
UDI004 DRR530881
UDI005 DRR530874
UDI006 DRR530876
UDI007 DRR530878
UDI008 DRR530880"

process() {
  bc=$1; run=$2
  line=$(awk -v r="$run" '$1==r' "$TMP/filereport.tsv")
  IFS=';' read -r url1 url2 <<< "$(cut -f2 <<< "$line")"
  IFS=';' read -r md51 md52 <<< "$(cut -f3 <<< "$line")"
  for r in 1 2; do
    url=url$r; md5=md5$r
    f="$TMP/${run}_${r}.fastq.gz"
    curl -sfL --retry 5 -C - -o "$f" "https://${!url}" || curl -sfL --retry 5 -o "$f" "https://${!url}"
    echo "${!md5}  $f" | md5sum -c --quiet
    seqtk sample -s"$SEED" "$f" "$N" | pigz -p 2 -c > "$OUT/SLX-11111.${bc}.AAA111.r_${r}.fq.gz.part"
    mv "$OUT/SLX-11111.${bc}.AAA111.r_${r}.fq.gz.part" "$OUT/SLX-11111.${bc}.AAA111.r_${r}.fq.gz"
    rm -f "$f"
  done
  echo "done $bc ($run)"
}
export -f process; export TMP OUT N SEED
printf "%s\n" "$MAP" | xargs -P 4 -L 1 bash -c 'process "$0" "$1"'
rm -rf "$TMP"
echo ALL_DONE
