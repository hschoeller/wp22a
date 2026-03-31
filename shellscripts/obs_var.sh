#!/bin/bash

#SBATCH --job-name=ObsVar
#SBATCH --output=./logs/ObsVar_%a.out
#SBATCH --error=./logs/ObsVar_%a.err
#SBATCH --array=1420
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=agpfahl
#SBATCH --mem=10G
#SBATCH --qos=agpfahl
#SBATCH --time=14-00:00:00

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_THREADS=1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

# cd /net/scratch/schoelleh96/WP2/WP2.2a
cd /home/schoelleh96/wp22a
source ~/miniforge3/etc/profile.d/conda.sh
conda activate wp22aR

# while true; do
#     used_ram=$(free -g | awk '/^Speicher:/ {print $3}')
#     echo "$(date): ${used_ram} GB" >> memory_usage.log
#     sleep 60
# done &

SEASON=("False") # ("False" "True")
WR=("True") #  "False")

for SEAS in "${SEASON[@]}"; do
  for WR_FLAG in "${WR[@]}"; do

    Rscript /home/schoelleh96/wp22a/RScripts/ObsVar.r \
        /scratch/schoelleh96/wp22a/ens_data/z_chunks/ ${SLURM_ARRAY_TASK_ID} \
        ${SEAS} ${WR_FLAG}

    # obsolete because ObsVar.r adds model data now
    # Rscript /home/schoelleh96/wp22a/RScripts/add_model_data.r \
    #     /scratch/schoelleh96/wp22a/ens_data/z_chunks/ ${SLURM_ARRAY_TASK_ID} \
    #     ${SEAS} ${WR_FLAG}

  done
done
