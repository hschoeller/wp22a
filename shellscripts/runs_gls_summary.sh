#!/bin/bash

#SBATCH --job-name=gls_summary
#SBATCH --output=./logs/summary_%A_%a.out
#SBATCH --error=./logs/summary_%A_%a.err
#SBATCH --partition=agpfahl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00
# SBATCH --array=0-3

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

Rscript /home/schoelleh96/wp22a/RScripts/run_gls_summary.r \
  /scratch/schoelleh96/wp22a/ens_data/z_chunks/ \
  "True" "False"
exit

SEASON=("False" "True")
WR=("True" "False")

# Expecting SLURM_ARRAY_TASK_ID in {0,1,2,3}
i=${SLURM_ARRAY_TASK_ID}

seas_idx=$(( i / 2 ))   # 0 or 1
wr_idx=$(( i % 2 ))     # 0 or 1

SEAS="${SEASON[$seas_idx]}"
WR_FLAG="${WR[$wr_idx]}"

Rscript /home/schoelleh96/wp22a/RScripts/run_gls_summary.r \
  /scratch/schoelleh96/wp22a/ens_data/z_chunks/ \
  "${SEAS}" "${WR_FLAG}"

