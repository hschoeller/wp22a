#!/bin/bash
#SBATCH --job-name=wcb_collect
#SBATCH --output=logs/wcb_collect_%j.out
#SBATCH --error=logs/wcb_collect_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1

CHUNK_DIR="/scratch/schoelleh96/wp22a/stage2_wcb_out_alt2"
OUT_FILE="/home/schoelleh96/wp22a/ens_data/wcb_plot_data_all_alt2.rds"


Rscript /home/schoelleh96/wp22a/RScripts/collect_wcb_plot_data_parallel.r \
  "$CHUNK_DIR" \
  "$OUT_FILE"