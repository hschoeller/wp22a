#!/bin/bash
#SBATCH --job-name=emmeans
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=logs/emmeans.out
#SBATCH --error=logs/emmeans.err
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

MODEL_DIR="/scratch/schoelleh96/wp22a/ens_data/z_chunks/modsSeasWR"
NCORES=$SLURM_CPUS_PER_TASK
OUT_FILE="/home/schoelleh96/wp22a/data/composites_emmeans.rds"

Rscript /home/schoelleh96/wp22a/RScripts/extract_wr_effects_emmeans.r "$MODEL_DIR" "$OUT_FILE" "$NCORES"