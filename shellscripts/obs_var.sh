#!/bin/bash

#SBATCH --job-name=ObsVar
#SBATCH --output=./ObsVar.log
#SBATCH --error=./ObsVar.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=main
#SBATCH --nodelist=calc03                                                                
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --mem=185G

# export R_LIBS=/net/scratch/schoelleh96/WP2/WP2.2a/RScripts/R_lib

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

Rscript /net/scratch/schoelleh96/WP2/WP2.2a/RScripts/ObsVar.r \
    /net/scratch/schoelleh96/WP2/WP2.2a/ens_data/zg_chunks/
