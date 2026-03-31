#!/bin/bash
#SBATCH --job-name=build_residual_cube
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=logs/build_residual_cube.out
#SBATCH --error=logs/build_residual_cube.err
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

MODEL_DIR="/scratch/schoelleh96/wp22a/ens_data/z_chunks/modsSeas"
CUBE="/home/schoelleh96/wp22a/ens_data/residual_cube_weighted.rds"
WEIGHT_MODE="stabilize"
NCORES=$SLURM_CPUS_PER_TASK
WRFILE="/home/schoelleh96/wp22a/data/wrnames.rds"
OUT_FILE="/home/schoelleh96/wp22a/data/composite_data.rds"

#Rscript /home/schoelleh96/wp22a/RScripts/build_residual_cube.r "$MODEL_DIR" "$CUBE" "$WEIGHT_MODE" "$NCORES"
Rscript /home/schoelleh96/wp22a/RScripts/compute_composites_from_cube.r "$CUBE" "$WRFILE" "$OUT_FILE" 1000 "$NCORES" 50