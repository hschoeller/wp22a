#!/bin/bash

#SBATCH --job-name=ObsVar
#SBATCH --output=./ObsVar.log
#SBATCH --error=./ObsVar.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=main
#SBATCH --nodelist=calc05                                                                   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --mem=375G

# export R_LIBS=/net/scratch/schoelleh96/WP2/WP2.2a/RScripts/R_lib

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