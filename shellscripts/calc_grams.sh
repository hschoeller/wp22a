#!/usr/bin/env bash

#SBATCH --job-name=calc_grams
#SBATCH --output=logs/grams_%j.out
#SBATCH --error=logs/grams_%j.err
#SBATCH --mem=250G
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

source ~/miniforge3/etc/profile.d/conda.sh
# conda activate wp22a

# SCRIPT="/home/schoelleh96/wp22a/pyscripts/calc_grams_anoms.py"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# python "${SCRIPT}"

export R_THREADS=1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

cd /home/schoelleh96/wp22a
conda activate wp22aR

Rscript RScripts/GramsAtt.r ens_data/residual_cube_weighted.rds \
    scratch/data/zg500_processed.nc data/wrnames.rds ./testGrams.rds \
    zg500_prime_lp 10000