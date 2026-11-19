#!/bin/bash
# Usage: sbatch submit_slurm.sh   (from the pipeline directory)
# The account is not stored here: set it once in ~/.bashrc with  export SBATCH_ACCOUNT=<your account>
# (or pass it at submission: sbatch -A <your account> submit_slurm.sh)

#SBATCH -J pipeline-bulkrnaseq-sm        # Job name
#SBATCH --mem=64G                         # Max memory (STAR needs ~60 GB per sample with a human index)
#SBATCH --time=2:00:00                    # Max runtime
#SBATCH --mail-type=ALL                   # Email notifications
#SBATCH -p sapphire                       # Partition
#SBATCH --cpus-per-task=16                # CPUs per task

# run.sh sizes STAR / FastQC from the CPUs and memory of this job
bash run.sh
