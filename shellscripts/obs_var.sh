#!/bin/bash

#SBATCH --job-name=ObsVar
#SBATCH --output=./ObsVar.log
#SBATCH --error=./ObsVar.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=main
#SBATCH --nodelist=calc04                                                                
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --mem=118G

# export R_LIBS=/net/scratch/schoelleh96/WP2/WP2.2a/RScripts/R_lib

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1

cd /net/scratch/schoelleh96/WP2/WP2.2a/RScripts

while true; do
    used_ram=$(free -g | awk '/^Speicher:/ {print $3}')
    echo "$(date): ${used_ram} GB" >> memory_usage.log
    sleep 60
done &

Rscript /net/scratch/schoelleh96/WP2/WP2.2a/RScripts/ObsVar.r

# source ~/mambaforge/bin/activate
# conda activate wp22a

# python ../pyscripts/splitter.py