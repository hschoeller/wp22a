#!/usr/bin/env bash

#SBATCH --job-name=upperwind_era5_allyears
#SBATCH --output=logs/wind_%j.out
#SBATCH --error=logs/wind_%j.err
#SBATCH --mem=200G
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00


source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22a

IN_DIR="/scratch/schoelleh96/wp22a/data/"
OUT_DIR="/scratch/schoelleh96/wp22a/data/upperwind"
SCRIPT="/home/schoelleh96/wp22a/pyscripts/compute_upperwind_dask.py"

mkdir -p "${OUT_DIR}" logs

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python "${SCRIPT}" \
  --in-dir "${IN_DIR}" \
  --out-dir "${OUT_DIR}" \
  --level 300 \
  --u-name "u_component_of_wind" \
  --v-name "v_component_of_wind" \
  --n-workers "${SLURM_CPUS_PER_TASK}" \
  --threads-per-worker 1 \
  --time-chunk 180 \
  --write-mode single