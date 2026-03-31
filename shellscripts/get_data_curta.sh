#!/bin/bash

#SBATCH --job-name=Get_Data
#SBATCH --output=logs/GetData_%a.out
#SBATCH --error=logs/GetData_%a.err

#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22a

# VARS=("geopotential" "geopotential" "total_cloud_cover" "medium_cloud_cover" "high_cloud_cover")
# VARS_SHORT=("z" "z" "tcc" "mcc" "hcc")
# DATASETS=("reanalysis-era5-pressure-levels"  "reanalysis-era5-pressure-levels" "reanalysis-era5-single-levels" "reanalysis-era5-single-levels" "reanalysis-era5-single-levels")
# PRODUCT_TYPES=("reanalysis" "ensemble_spread" "reanalysis" "reanalysis" "reanalysis")
# ENS_VALUES=("False" "True" "False" "False" "False")

# i=$SLURM_ARRAY_TASK_ID
# VAR="${VARS[$i]}"
# VAR_SHORT="${VARS_SHORT[$i]}"
# DATASET="${DATASETS[$i]}"
# PRODUCT_TYPE="${PRODUCT_TYPES[$i]}"
# ENS="${ENS_VALUES[$i]}"

DATASET="reanalysis-era5-pressure-levels"
PRODUCT_TYPE="reanalysis"
ENS="False"
OUTROOT="/scratch/schoelleh96/wp22a"

# TASKS=(
#   "u  u_component_of_wind 500"
#   "v  v_component_of_wind 500"
#   "t  temperature         500"
#   "u  u_component_of_wind 850"
#   "v  v_component_of_wind 850"
#   "t  temperature         850"
#   "z  geopotential        850"
#   "u  u_component_of_wind 300"
#   "v  v_component_of_wind 300"
# )

# TASKS=(
#   "r  relative_humidity 850"
#   "r  relative_humidity 825"
#   "r  relative_humidity 800"
#   "r  relative_humidity 775"
#   "r  relative_humidity 750"
#   "r  relative_humidity 700"
# )

# read -r VAR_SHORT VAR PLEV <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"

DATASET="reanalysis-era5-single-levels"
ENS="False"
PRODUCT_TYPE="reanalysis"

TASKS=(
  "lsm  land_sea_mask"
  "tclw  total_column_cloud_liquid_water"
  "tciw  total_column_cloud_ice_water"
  "tcrw  total_column_rain_water"
  "tcsw  total_column_snow_water"
)

read -r VAR_SHORT VAR  <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"


if [ "$ENS" = "True" ]; then
    BASEFOLDER="ens_data"
else
    BASEFOLDER="data"
fi

python ../pyscripts/get_data.py "${VAR}" "${DATASET}" "${PRODUCT_TYPE}" "${ENS}" \
    "${VAR_SHORT}" "/scratch/schoelleh96/wp22a" # "${PLEV}"
# python ../pyscripts/splitter.py "/scratch/schoelleh96/wp22a/${BASEFOLDER}/${VAR}.nc" \
#     "/scratch/schoelleh96/wp22a/${BASEFOLDER}/${VAR_SHORT}_chunks/" --chunks 11 241


