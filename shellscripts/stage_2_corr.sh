#!/bin/bash

#SBATCH --job-name=Stage2Corr
#SBATCH --output=./logs/Stage2Corr_%a.out
#SBATCH --error=./logs/Stage2Corr_%a.err
#SBATCH --array=0-2652%500
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=3G
#SBATCH --qos=agpfahl
#SBATCH --time=2-00:00:00

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

CHUNK_DIR="/scratch/schoelleh96/wp22a/stage2_chunks"
OUT_BASE="/scratch/schoelleh96/wp22a/stage2_correlations"
SCRIPT="/home/schoelleh96/wp22a/RScripts/stage2_corr.r"

mkdir -p ./logs
mkdir -p "${OUT_BASE}"

Rscript "${SCRIPT}" \
    "${CHUNK_DIR}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${OUT_BASE}" \
    "true"