#!/bin/bash
#SBATCH --job-name=psd_mean_iqr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=logs/psd_mean_iqr_%j.out
#SBATCH --error=logs/psd_mean_iqr_%j.err
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl

set -euo pipefail
mkdir -p logs

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

# Avoid oversubscription (FFT/BLAS)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1


CHUNK_DIR="/scratch/schoelleh96/wp22a/ens_data/z_chunks"
NC_VAR="z"
OUT_DIR="/home/schoelleh96/wp22a/ens_data/psd_out"

# Use all allocated cores, but script internally caps chunk cores to 32 by default
Rscript /home/schoelleh96/wp22a/RScripts/CalcPSD.r "${CHUNK_DIR}" "${NC_VAR}" "${OUT_DIR}" "${SLURM_CPUS_PER_TASK}"
