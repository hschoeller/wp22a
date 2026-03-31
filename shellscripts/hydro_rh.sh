#!/bin/bash
#SBATCH --job-name=hydro_rh
#SBATCH --output=logs/hydro_rh_%j.out
#SBATCH --error=logs/hydro_rh_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22a

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# hydrometeor sum
python /home/schoelleh96/wp22a/pyscripts/hydrosum.py \
    /scratch/schoelleh96/wp22a/data/total_column_snow_waterNone.nc \
    /scratch/schoelleh96/wp22a/data/total_column_rain_waterNone.nc \
    /scratch/schoelleh96/wp22a/data/total_column_cloud_liquid_waterNone.nc \
    /scratch/schoelleh96/wp22a/data/total_column_cloud_ice_waterNone.nc

# RH layer mean
python /home/schoelleh96/wp22a/pyscripts/rh_layer_mean.py \
    /scratch/schoelleh96/wp22a/data/relative_humidity700.nc \
    /scratch/schoelleh96/wp22a/data/relative_humidity750.nc \
    /scratch/schoelleh96/wp22a/data/relative_humidity775.nc \
    /scratch/schoelleh96/wp22a/data/relative_humidity800.nc \
    /scratch/schoelleh96/wp22a/data/relative_humidity825.nc \
    /scratch/schoelleh96/wp22a/data/relative_humidity850.nc

