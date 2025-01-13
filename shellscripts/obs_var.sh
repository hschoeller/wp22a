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
#SBATCH --mem=180G

export R_LIBS=/net/scratch/schoelleh96/WP2/WP2.2a/RScripts/R_lib

cd /net/scratch/schoelleh96/WP2/WP2.2a/RScripts

Rscript /net/scratch/schoelleh96/WP2/WP2.2a/RScripts/ObsVar.r

