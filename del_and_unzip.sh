#!/bin/bash
#SBATCH --job-name=unzip_job
#SBATCH --output=unzip_%j.out
#SBATCH --error=unzip_%j.err
#SBATCH --partition=agpfahl
#SBATCH --qos=agpfahl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=02:00:00

set -euo pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

### ========== Edit these three paths ==========
ZIP="/scratch/schoelleh96/wp22a/ens_data/zg_mods.zip"
TARGET="/scratch/schoelleh96/wp22a/ens_data/"    # where files ultimately should live
REMOVE1="/scratch/schoelleh96/wp22a/ens_data/mods"
REMOVE2="/scratch/schoelleh96/wp22a/ens_data/modsbackup"
### =============================================

echo "Job started on $(hostname) at $(date)"
echo "ZIP: $ZIP"
echo "TARGET: $TARGET"
echo "Removing: $REMOVE1 $REMOVE2"

# 1) Remove the two target folders (be careful!)
rm -rf "$REMOVE1" "$REMOVE2"
echo "Removed old folders."

unzip -q "$ZIP" -d "$TARGET"
echo "Extraction finished in $LOCAL_OUT"

echo "Completed at $(date)"
