#!/bin/bash

#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl
#SBATCH --job-name=z500_diag
#SBATCH --output=logs/z500_diag_%j.out
#SBATCH --error=logs/z500_diag_%j.err

mkdir -p logs

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22a

# Optional: thread settings
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- Paths ---
PYTHON_SCRIPT=/home/schoelleh96/wp22a/pyscripts/calc_grad_lap.py
INPUT_FILE=/scratch/schoelleh96/wp22a/data/geopotential500_full.nc
GRAD_OUTPUT_FILE=/scratch/schoelleh96/wp22a/data/z500_grad_mag.nc
LAP_OUTPUT_FILE=/scratch/schoelleh96/wp22a/data/z500_laplacian.nc
ABS_LAP_OUTPUT_FILE=/scratch/schoelleh96/wp22a/data/z500_abs_laplacian.nc

# --- Run ---
python "${PYTHON_SCRIPT}" \
    --input "${INPUT_FILE}" \
    --grad-output "${GRAD_OUTPUT_FILE}" \
    --lap-output "${LAP_OUTPUT_FILE}" \
    --abs-lap-output "${ABS_LAP_OUTPUT_FILE}" \
    --var z 
    # --convert-geopotential-to-height \
    # --lat latitude \
    # --lon longitude \
    # --time time \
    # --chunks time=20
