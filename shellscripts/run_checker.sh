#!/bin/bash

#SBATCH --job-name=checkerboard
#SBATCH --output=logs/checkerboard.out
#SBATCH --error=logs/checkerboard.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=agpfahl
#SBATCH --mail-type=ALL
#SBATCH --mem=30G
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22a

INFILE=/scratch/schoelleh96/wp22a/ens_data/log_geopotential.nc
OUTDIR=/home/schoelleh96/wp22a/checker_detect_log
mkdir -p "${OUTDIR}"

python3 /home/schoelleh96/wp22a/pyscripts/detect_checkerboard.py "${INFILE}" "${OUTDIR}" --plot_per_lat