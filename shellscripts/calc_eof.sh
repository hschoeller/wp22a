#!/bin/bash

#SBATCH --job-name=eofs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --nodelist=calc05
#SBATCH --output=./eof.out                                                                    
#SBATCH --error=./eof.err                                                                     
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de    
#SBATCH --array=1-9

# This script will calculate EOFs based on the gh500 anomalies in the Euro-Atlantic region

echo "SLURM CLUSTER NAME: $SLURM_CLUSTER_NAME"
echo "SLURM CPUS PER TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM JOB ACCOUNT : $SLURM_JOB_ACCOUNT"
echo "SLURM JOB ID: $SLURM_JOB_ID"
echo "slurm job name: $SLURM_JOB_NAME"
echo "slurm job nodelist: $SLURM_JOB_NODELIST"
echo "slurm job partition: $SLURM_JOB_PARTITION"
echo "slurm job uid: $SLURM_JOB_UID"
echo "slurm job user: $SLURM_JOB_USER"
echo "slurm job procid : $SLURM_PROCID"
echo "slurm step id: $SLURM_STEP_ID"
echo "slurm step num tasks: $SLURM_STEP_NUM_TASKS"
echo "Running on $(hostname) with $SLURM_CPUS_PER_TASK CPUs"

# Modify OMP_NUM_THREADS to reflect available CPUs
# export OMP_NUM_THREADS=48

echo "Calculating EOFs"

source ~/mambaforge/bin/activate
conda activate wp22a

python ../pyscripts/prepare_eof.py $SLURM_ARRAY_TASK_ID

# conda activate eof_env
# python ../pyscripts/do_eof.py

#export CDO_WEIGHT_MODE=on
#cdo -O -P 48 -eoftime,50 ../data/zg_norm.nc ../data/evals.nc ../data/evecs.nc
#cdo -O -P 48 -eofcoeff ../data/evecs.nc ../data/zg_norm.nc ../data/pc
#cdo -O -merge ../data/pc0* ../data/pcs.nc
#rm -f ../data/pc0*.nc



