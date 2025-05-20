#!/bin/bash

#SBATCH --job-name=CompPermute
#SBATCH --output=./CompPermute.log
#SBATCH --error=./CompPermute.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=main
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --mem=118G
#SBATCH --nodelist=calc04                                                                  

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1

cd /net/scratch/schoelleh96/WP2/WP2.2a

while true; do
    used_ram=$(free -g | awk '/^Speicher:/ {print $3}')
    echo "$(date): ${used_ram} GB" >> memory_usage.log
    sleep 60
done &

# Rscript /net/scratch/schoelleh96/WP2/WP2.2a/RScripts/composites_permuted.r
Rscript /net/scratch/schoelleh96/WP2/WP2.2a/RScripts/WR_Log_Var_Raw.r