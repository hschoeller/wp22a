#!/bin/bash

#SBATCH --job-name=Get_Data
#SBATCH --output=logs/GetData_%A_%a.out
#SBATCH --error=logs/GetData_%A_%a.err

#SBATCH --array=1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de
#SBATCH --mem=40G

source ~/mambaforge/bin/activate
conda activate wp22a

VARS=("geopotential" "geopotential" "total_cloud_cover" "medium_cloud_cover" "high_cloud_cover")
VARS_SHORT=("zg" "zg" "tcc" "mcc" "hcc")
DATASETS=("reanalysis-era5-pressure-levels"  "reanalysis-era5-pressure-levels" "reanalysis-era5-single-levels" "reanalysis-era5-single-levels" "reanalysis-era5-single-levels")
PRODUCT_TYPES=("reanalysis" "ensemble_spread" "reanalysis" "reanalysis" "reanalysis")
ENS_VALUES=("False" "True" "False" "False" "False")

i=$SLURM_ARRAY_TASK_ID
VAR="${VARS[$i]}"
VAR_SHORT="${VARS_SHORT[$i]}"
DATASET="${DATASETS[$i]}"
PRODUCT_TYPE="${PRODUCT_TYPES[$i]}"
ENS="${ENS_VALUES[$i]}"

if [ "$ENS" = "True" ]; then
    BASEFOLDER="ens_data"
else
    BASEFOLDER="data"
fi

python ../pyscripts/get_data.py "${VAR}" "${DATASET}" "${PRODUCT_TYPE}" "${ENS}" \
    "/net/scratch/schoelleh96/WP2/WP2.2a/"
# python ../pyscripts/splitter.py "/net/scratch/schoelleh96/WP2/WP2.2a/${BASEFOLDER}/${VAR}.nc" \
#     "/net/scratch/schoelleh96/WP2/WP2.2a/${BASEFOLDER}/${VAR_SHORT}_chunks/" --chunks 1 48 --no-regrid
