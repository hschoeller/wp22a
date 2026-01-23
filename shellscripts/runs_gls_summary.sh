#!/bin/bash

#SBATCH --job-name=gls_summary
#SBATCH --output=./logs/summary.out
#SBATCH --error=./logs/summary.err
#SBATCH --partition=agpfahl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

Rscript /home/schoelleh96/wp22a/RScripts/run_gls_summary.r
