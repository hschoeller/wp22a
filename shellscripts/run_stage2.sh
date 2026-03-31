#!/bin/bash

#SBATCH --job-name=Stage2
#SBATCH --output=./logs/Stage2_%a.out
#SBATCH --error=./logs/Stage2_%a.err
#SBATCH --array=0-10%10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=12G
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

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
OUT_BASE="/scratch/schoelleh96/wp22a/stage2_model_estimates"
WR_RDS="/home/schoelleh96/wp22a/data/wrnames.rds"
SCRIPT="/home/schoelleh96/wp22a/RScripts/stage2Mods.r"

mkdir -p ./logs
mkdir -p "${OUT_BASE}"

# first pass: simple single-predictor models for plots 1/2
Rscript "${SCRIPT}" \
    "${CHUNK_DIR}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${OUT_BASE}" \
    "${WR_RDS}" \
    "false"


# second pass: only run this after you have decided which predictors to keep for plots 3/4
# Rscript "${SCRIPT}" \
#     "${CHUNK_DIR}" \
#     "${OUT_BASE}/full" \
#     "${WR_RDS}" \
#     "full"