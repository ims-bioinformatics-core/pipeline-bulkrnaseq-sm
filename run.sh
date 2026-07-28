#!/bin/bash

# ==========================
# Pipeline directories
# ==========================
PIPELINE_DIR="/rds/project/rds-O11U8YqSuCk/bioinformatics/pipelines/pipeline-bulkrnaseq-sm"
CONFIG_FILE="$PIPELINE_DIR/config.yml"
# ==========================
# Run Snakemake
# ==========================
start=$(date +%s)
snakemake -s "$PIPELINE_DIR/Snakefile" \
          -j 2 \
          --use-conda \
          --rerun-incomplete \
          --configfile "$CONFIG_FILE" \
          --config star_threads=4 fastqc_threads=4 pipeline_dir="$PIPELINE_DIR"

end=$(date +%s)
runtime=$((end - start))
runtime_minutes=$(echo "scale=2; $runtime/60" | bc)
echo "Pipeline finished in $runtime seconds ($runtime_minutes minutes)"
