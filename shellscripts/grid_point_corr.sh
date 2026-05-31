#!/bin/bash

#SBATCH --job-name=Grid_Corr
#SBATCH --output=./logs/Grid_Corr_%a.out
#SBATCH --error=./logs/Grid_Corr_%a.err
#SBATCH --array=0 # 8-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=100G
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR


RESID_FILE="/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
DATA_DIR="/scratch/schoelleh96/wp22a/data"
OUT_DIR="/home/schoelleh96/wp22a/ens_data/corr_maps"


COV_FILES=(
  "$DATA_DIR/eady.nc"
  "$DATA_DIR/medium_cloud_cover.nc"
  "$DATA_DIR/high_cloud_cover.nc"
  "$DATA_DIR/hydrosum.nc"
  "$DATA_DIR/relative_humidity_500-650.nc"
  "$DATA_DIR/relative_humidity_500-850.nc"
  "$DATA_DIR/relative_humidity_700-850.nc"
  "$DATA_DIR/wind.nc"
  "$DATA_DIR/z500_grad_mag.nc"
  "$DATA_DIR/z500_laplacian.nc"
  "$DATA_DIR/z500_abs_laplacian.nc"
)

COV_VARS=(
  "eady_growth_rate"
  "mcc"
  "hcc"
  "hydrosum"
  "r"
  "r"
  "r"
  "upper_wind_speed"
  "z_grad_mag"
  "z_laplacian"
  "z_abs_laplacian"
)

START_YEARS=(1950 1979 1950)
END_YEARS=(1979 2024 2024)

i=$SLURM_ARRAY_TASK_ID

COV_FILE="${COV_FILES[$i]}"
COV_VAR="${COV_VARS[$i]}"

COV_FILE="$DATA_DIR/zg500_processed.nc"
COV_VAR="zg500_prime_lp"

for j in "${!START_YEARS[@]}"; do
  START="${START_YEARS[$j]}"
  END="${END_YEARS[$j]}"
  COV_NAME=$(basename "$COV_FILE" .nc)
OUT_FILE="$OUT_DIR/${COV_NAME}_${START}_${END}_corr_map.rds"
  Rscript /home/schoelleh96/wp22a/RScripts/corr_var_gridpoints.r \
    "$RESID_FILE" \
    "$COV_FILE"   \
    "$COV_VAR"    \
    "$OUT_FILE"   \
    "$START"      \
    "$END"
done
