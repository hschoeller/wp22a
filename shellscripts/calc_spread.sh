#!/bin/bash

#SBATCH --job-name=ensemble_spread
#SBATCH --output=./spread_calc.log
#SBATCH --error=./spread_calc.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --nodelist=calc04                                                                  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --array=82-83

source ~/mambaforge/bin/activate
conda activate wp22a

python ../pyscripts/calc_spread.py $SLURM_ARRAY_TASK_ID
