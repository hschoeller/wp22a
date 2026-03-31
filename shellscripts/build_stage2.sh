#!/usr/bin/env bash
#SBATCH --job-name=build_stage2
#SBATCH --output=logs/build_stage2_%a.out
#SBATCH --error=logs/build_stage2_%a.err
#SBATCH --array=1513
#SBATCH --time=14-00:00:00
#SBATCH --mem=9G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl


cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1

RESP_RDS="/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
MCC_NC="/scratch/schoelleh96/wp22a/data/medium_cloud_cover.nc"
HCC_NC="/scratch/schoelleh96/wp22a/data/high_cloud_cover.nc"
WCB_ROOT="/scratch/schoelleh96/wp22a/ELIAS_data"
EADY_NC="/scratch/schoelleh96/wp22a/data/eady.nc"
WIND_NC="/scratch/schoelleh96/wp22a/data/wind.nc"
HYDRO_NC="/scratch/schoelleh96/wp22a/data/hydrosum.nc"
RH_NC="/scratch/schoelleh96/wp22a/data/relative_humidity_700-850.nc"
GRAD_NC=/scratch/schoelleh96/wp22a/data/z500_grad_mag.nc
LAP_NC=/scratch/schoelleh96/wp22a/data/z500_laplacian.nc

MCC_VAR="mcc"
HCC_VAR="hcc"
EADY_VAR="eady_growth_rate"
WIND_VAR="upper_wind_speed"
HYDRO_VAR="hydrosum"
RH_VAR="r"
GRAD_VAR="z_grad_mag"
LAP_VAR="z_laplacian"

R_SCRIPT="/home/schoelleh96/wp22a/RScripts/build_stage2_chunks.r"
OUT_DIR="/scratch/schoelleh96/wp22a/stage2_chunks"

LAT_CHUNKS=11
LON_CHUNKS=241
OVERWRITE=0

mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/logs"

Rscript "${R_SCRIPT}" \
  --task-id="${SLURM_ARRAY_TASK_ID}" \
  --chunks="${LAT_CHUNKS},${LON_CHUNKS}" \
  --response-rds="${RESP_RDS}" \
  --mcc-nc="${MCC_NC}" \
  --hcc-nc="${HCC_NC}" \
  --mcc-var="${MCC_VAR}" \
  --hcc-var="${HCC_VAR}" \
  --eady-nc="${EADY_NC}" \
  --wind-nc="${WIND_NC}" \
  --eady-var="${EADY_VAR}" \
  --wind-var="${WIND_VAR}" \
  --rh-nc="${RH_NC}" \
  --hydo_nc="${HYDRO_NC}" \
  --rh_var="${RH_VAR}" \
  --hydro_var="${HYDRO_VAR}" \
  --grad-nc="${GRAD_NC}" \
  --lap_nc="${LAP_NC}" \
  --grad_var="${GRAD_VAR}" \
  --lap_var="${LAP_VAR}" \
  --wcb-root="${WCB_ROOT}" \
  --out-dir="${OUT_DIR}" \
  --overwrite="${OVERWRITE}"