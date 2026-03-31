#!/bin/bash

#SBATCH --job-name=test
#SBATCH --output=./logs/test%a.out
# SBATCH --error=./logs/test%a.err
# SBATCH --array=0-2650
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mail-type=ALL
#SBATCH --mem=2G
#SBATCH --qos=standard
#SBATCH --time=14-00:00:00

cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

Rscript /home/schoelleh96/wp22a/RScripts/test.r
