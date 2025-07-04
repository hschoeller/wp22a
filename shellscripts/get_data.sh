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

VARS=("geopotential" "total_cloud_cover" "medium_cloud_cover" "high_cloud_cover")
VARS_SHORT=("zg" "tcc" "mcc" "hcc")
ENS_VALUES=("True" "False" "False" "False")

for i in "${!VARS[@]}"; do
    VAR="${VARS[$i]}"
    VAR_SHORT="${VARS_SHORT[$i]}"
    ENS="${ENS_VALUES[$i]}"

    if [ "$ENS" = "True" ]; then
        BASEFOLDER="ens_data"
    else
        BASEFOLDER="data"
    fi

    # python ../pyscripts/get_data.py "${VAR}" "${ENS}"
    python ../pyscripts/splitter.py "../${BASEFOLDER}/${VAR}.nc" \
        "../${BASEFOLDER}/${VAR_SHORT}_chunks/" --chunks 1 48 --no-regrid
done