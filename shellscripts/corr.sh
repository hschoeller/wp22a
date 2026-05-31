#!/bin/bash

#SBATCH --job-name=yearcorr
#SBATCH --output=./logs/yearcorr_%A_%a.out
#SBATCH --error=./logs/yearcorr_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=60G
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00
#SBATCH --array=1

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

SCRIPT="/home/schoelleh96/wp22a/RScripts/yearly_grid_corrs.r"
RESP_RDS="/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
OUT_DIR="/home/schoelleh96/wp22a/ens_data/yearly_corrs"

mkdir -p "${OUT_DIR}" ./logs

# NAMES=(
#   "mcc"
#   "hcc"
#   "eady"
#   "wind"
#   "hydrosum"
  # "rh"
# )

# NC_FILES=(
#   "/scratch/schoelleh96/wp22a/data/medium_cloud_cover.nc"
#   "/scratch/schoelleh96/wp22a/data/high_cloud_cover.nc"
#   "/scratch/schoelleh96/wp22a/data/eady.nc"
#   "/scratch/schoelleh96/wp22a/data/wind.nc"
#   "/scratch/schoelleh96/wp22a/data/hydrosum.nc"
  # "/scratch/schoelleh96/wp22a/data/relative_humidity_500-850.nc"
# )

# VAR_NAMES=(
#   "mcc"
#   "hcc"
#   "eady_growth_rate"
#   "upper_wind_speed"
#   "hydrosum"
#   "r"
# )

NAMES=(
  "z_grad_mag"
  "z_laplacian"
  "z_abs_laplacian"
)

NC_FILES=(
  "/scratch/schoelleh96/wp22a/data/z500_grad_mag.nc"
  "/scratch/schoelleh96/wp22a/data/z500_laplacian.nc"
  "/scratch/schoelleh96/wp22a/data/z500_abs_laplacian.nc"
)

VAR_NAMES=(
  "z_grad_mag"
  "z_laplacian"
  "z_abs_laplacian"
)


i=${SLURM_ARRAY_TASK_ID}

NAME="${NAMES[$i]}"
NC_FILE="${NC_FILES[$i]}"
VAR_NAME="${VAR_NAMES[$i]}"
OUT_FILE="${OUT_DIR}/yearly_corr_${NAME}.rds"

echo "Running variable: ${NAME}"
echo "NC file: ${NC_FILE}"
echo "Variable: ${VAR_NAME}"
echo "Output: ${OUT_FILE}"

NC_FILE="/scratch/schoelleh96/wp22a/data/zg500_processed.nc"
VAR_NAME="zg500_prime_lp"
OUT_FILE="${OUT_DIR}/yearly_corr_zg500_lp.rds"

Rscript "${SCRIPT}" \
  "${RESP_RDS}" \
  "${NC_FILE}" \
  "${VAR_NAME}" \
  "${OUT_FILE}"