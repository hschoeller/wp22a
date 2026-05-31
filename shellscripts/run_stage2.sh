#!/bin/bash

#SBATCH --job-name=Stage2
#SBATCH --output=./logs/Stage2_%a.out
#SBATCH --error=./logs/Stage2_%a.err
# SBATCH --array=1-2651
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=2G
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

# CHUNK_DIR="/scratch/schoelleh96/wp22a/stage2_chunks"
# OUT_BASE="/scratch/schoelleh96/wp22a/stage2_model_estimates"
# WR_RDS="/home/schoelleh96/wp22a/data/wrnames.rds"
# SCRIPT="/home/schoelleh96/wp22a/RScripts/stage2Mods.r"

CHUNK_DIR="/scratch/schoelleh96/wp22a/stage2_chunks_small"
# OUT_BASE="/scratch/schoelleh96/wp22a/stage2_models_small"
OUT_BASE="/scratch/schoelleh96/wp22a/stage1_models_perm"
WR_RDS="/home/schoelleh96/wp22a/data/wrnames.rds"
SCRIPT="/home/schoelleh96/wp22a/RScripts/stage1_perm.r"

mkdir -p ./logs
mkdir -p "${OUT_BASE}"

Rscript "${SCRIPT}" \
    "${CHUNK_DIR}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${OUT_BASE}" \
    "${WR_RDS}" \
    "true" \
    "10000"


