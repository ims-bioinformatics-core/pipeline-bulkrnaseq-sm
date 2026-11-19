#!/bin/bash
# Defaults for `sbatch run.sh`; options on the sbatch command line override them.
# The account is not stored here: set it once in ~/.bashrc with  export SBATCH_ACCOUNT=<your account>
# STAR needs ~60 GB per sample (human index), so the job needs at least 64 GB.
#SBATCH -J pipeline-bulkrnaseq-sm
#SBATCH -p icelake
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00

# ==========================
# Pipeline directories
# ==========================
PIPELINE_DIR="/rds/project/rds-O11U8YqSuCk/bioinformatics/pipelines/pipeline-bulkrnaseq-sm"
CONFIG_FILE="$PIPELINE_DIR/config.yml"

# ==========================
# Resources
# ==========================
# SLURM_CPUS_ON_NODE is what the job really has: asking for more memory makes SLURM add CPUs.
CORES="${SLURM_CPUS_ON_NODE:-$(nproc)}"
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    TOTAL_MEM_MB="$SLURM_MEM_PER_NODE"
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    # Job submitted with --mem-per-cpu
    TOTAL_MEM_MB=$((SLURM_MEM_PER_CPU * CORES))
elif [ -n "${SLURM_JOB_ID:-}" ]; then
    # Job submitted without any memory option: SLURM sets no variable, ask it for the job's limit
    TOTAL_MEM_MB=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | awk -v cores="$CORES" '
        function mb(v) { u = substr(v, length(v)); n = v + 0;
                         return u == "T" ? n * 1048576 : u == "G" ? n * 1024 : u == "K" ? n / 1024 : n }
        { for (i = 1; i <= NF; i++) {
              if ($i ~ /^MinMemoryNode=/) { split($i, a, "="); print int(mb(a[2])); exit }
              if ($i ~ /^MinMemoryCPU=/)  { split($i, a, "="); print int(mb(a[2]) * cores); exit } } }')
    TOTAL_MEM_MB="${TOTAL_MEM_MB:-4096}"
else
    TOTAL_MEM_MB=$(awk '/MemAvailable/ {print int($2 / 1024)}' /proc/meminfo)
fi
MEM_MB=$((TOTAL_MEM_MB * 97 / 100))

# One STAR job peaks at ~58 GB with a human index (genome in memory + BAM sorting), so memory
# decides how many run at once; the CPUs are then shared between them.
STAR_MEM=60                  # GB per STAR job
STAR_MAX_THREADS=24          # STAR gains little beyond this
LIGHT_JOB_MEM_MB=16000       # every other rule: limits how many run at once

STAR_JOBS=$((MEM_MB / (STAR_MEM * 1024)))
[ "$STAR_JOBS" -lt 1 ] && STAR_JOBS=1
STAR_THREADS=$((CORES / STAR_JOBS))
[ "$STAR_THREADS" -gt "$STAR_MAX_THREADS" ] && STAR_THREADS="$STAR_MAX_THREADS"
[ "$STAR_THREADS" -lt 1 ] && STAR_THREADS=1

echo "Resources: $CORES cores, $MEM_MB MB -> up to $STAR_JOBS STAR samples at once, $STAR_THREADS threads / $STAR_MEM GB each"
if [ "$MEM_MB" -lt $((STAR_MEM * 1024)) ]; then
    echo "WARNING: $MEM_MB MB is less than one STAR job needs ($STAR_MEM GB): Snakemake will not start STAR."
    echo "         Submit with more memory, e.g. sbatch -p icelake --cpus-per-task=16 --mem=64G run.sh"
fi

# ==========================
# Run Snakemake
# ==========================
start=$(date +%s)
snakemake -s "$PIPELINE_DIR/Snakefile" \
          --cores "$CORES" \
          --resources mem_mb="$MEM_MB" \
          --default-resources mem_mb="$LIGHT_JOB_MEM_MB" \
          --use-conda \
          --rerun-incomplete \
          --configfile "$CONFIG_FILE" \
          --config star_threads="$STAR_THREADS" star_mem="$STAR_MEM" fastqc_threads=4 pipeline_dir="$PIPELINE_DIR" \
          "$@"

end=$(date +%s)
runtime=$((end - start))
runtime_minutes=$(echo "scale=2; $runtime/60" | bc)
echo "Pipeline finished in $runtime seconds ($runtime_minutes minutes)"
