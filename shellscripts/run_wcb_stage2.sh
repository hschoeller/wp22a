#!/bin/bash
#SBATCH --job-name=wcb_stage2
#SBATCH --output=./logs/wcb_stage2_%a.out
#SBATCH --error=./logs/wcb_stage2_%a.err
#SBATCH --array=1-2651%500
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=6G
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

CHUNK_DIR="/scratch/schoelleh96/wp22a/stage2_chunks_res"
OUT_DIR="/scratch/schoelleh96/wp22a/stage2_wcb_out_alt2"
WCB_RDS="/home/schoelleh96/wp22a/data/wcb_predictors_all_years.rds"

# # Run the R script for this array index
Rscript /home/schoelleh96/wp22a/RScripts/wcb_stage2.r \
  "$CHUNK_DIR" \
  "$SLURM_ARRAY_TASK_ID" \
  "$OUT_DIR" \
  "$WCB_RDS" \
  "inf_12,asc_00,out_12" \
  true


# Rscript /home/schoelleh96/wp22a/RScripts/wcb_wr_mods.r \
#   --task-id=${SLURM_ARRAY_TASK_ID} \
#   --in-dir=/scratch/schoelleh96/wp22a/stage2_chunks_wcb \
#   --out-dir=/scratch/schoelleh96/wp22a/wcb_mod_local \
#   --wr-rds=/home/schoelleh96/wp22a/data/wrnames.rds \
#   --use-wr=0

# Rscript /home/schoelleh96/wp22a/RScripts/wr_wcb_stacked.r \
#   --task-id=${SLURM_ARRAY_TASK_ID} \
#   --in-dir=/scratch/schoelleh96/wp22a/stage2_chunks_wcb \
#   --out-dir=/scratch/schoelleh96/wp22a/wr_wcb_stacked \
#   --wr-rds=/home/schoelleh96/wp22a/data/wrnames.rds