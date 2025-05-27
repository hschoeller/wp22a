#!/bin/bash

#SBATCH --job-name=Get_Data
#SBATCH --output=./GetData.log
#SBATCH --error=./GetData.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=main
#SBATCH --nodelist=calc03                                                             
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --mem=180G

source ~/mambaforge/bin/activate
conda activate wp22a

# python ../pyscripts/get_data.py
python ../pyscripts/splitter.py ../data/total_cloud_cover.nc \
    ../data/tcc_chunks/ --chunks 1 48 --regrid