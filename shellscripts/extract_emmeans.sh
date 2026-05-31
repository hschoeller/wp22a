#!/bin/bash
#SBATCH --job-name=extract_emm_chunks
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=logs/extract_emm_chunks%j.out
#SBATCH --error=logs/extract_emm_chunks%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

STAGE2_DIR="/scratch/schoelleh96/wp22a/stage2_model_estimates"
OUTFILE="/home/schoelleh96/wp22a/ens_data/wr_2ndemmeans.rds"

Rscript ./RScripts/extract_wr_emmeans_from_chunks.r \
    "$STAGE2_DIR" \
    "$OUTFILE" \
    "$SLURM_CPUS_PER_TASK"