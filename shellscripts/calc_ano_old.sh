#!/bin/bash

#SBATCH --job-name=anos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=main
#SBATCH --nodelist=calc04
#SBATCH --output=./anos.out
#SBATCH --error=./anos.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schoeller@fu-berlin.de

# This script will calculate daily gh500 means for each era5 grid point inside the Euro-Atlantic region

#source /home/schoelleh96/mambaforge/condabin/mamba

#mamba activate wp22a

#<<'COMMENT'

ERA5_DIR=/daten/reana/arch/reanalysis/reanalysis/DKRZ/IFS/ERA5/day/atmos/zg/r1i1p1-050000Pa/

echo "Calculating mean"

N=1
MERGE=''
for FILE in ${ERA5_DIR}zg_day_reanalysis_era5_r1i1p1-050000Pa_{1950..2023}*.nc
do
    echo $FILE
    NEW_FILE=subset_$N
    cdo -sellonlatbox,-110,0,20,90 ${FILE} ${NEW_FILE} && \
        MERGE="${MERGE} ${NEW_FILE}"
    N=$((N+1))
done
cdo -mergetime ${MERGE} temporally_merged.nc
cdo -ydaymean temporally_merged.nc ../data/ydailymean.nc

echo "Calulating anomalies"

N=1
for Y in {1950..2023}
do
    echo "Year=$Y"
    NEW_FILE=subset_$N
    cdo -sub ${NEW_FILE} ../data/ydailymean.nc ../data/${Y}_ano.nc
    N=$((N+1))
done

rm -f ${MERGE}
rm -f temporally_merged.nc

cdo -O -mergetime ../data/*_ano.nc ../data/anos.nc

rm -f ../data/*_ano.nc

cdo -P 48 -detrend ../data/anos.nc ../data/anosdetrend.nc

echo "Calculating moving window std"

cdo runstd,30 ../data/anosdetrend.nc ../data/std30.nc
cdo settime,12:00 ../data/std30.nc ../data/std30_12h.nc
cdo fldmean ../data/std30_12h.nc ../data/std30fldmean.nc

echo "Removing seasonal cycle by normalization with std"

cdo seldate,1950-01-16,2023-12-17 ../data/anosdetrend.nc ../data/anosdetrend_trimmed.nc

# COMMENT

###

cdo div ../data/anosdetrend_trimmed.nc -enlarge,../data/anosdetrend_trimmed.nc ../data/std30fldmean.nc ../data/anos_norm.nc

#cdo griddes ../data/anosdetrend_trimmed.nc > ../data/grid.txt

#cdo remapbil,../data/grid.txt ../data/std30fldmean.nc ../data/std30fldmean_remap.nc

#cdo div ../data/anosdetrend_trimmed.nc ../data/std30fldmean_remap.nc ../data/anos_norm.nc

#python ../pyscripts/normalize.py ../data/anosdetrend_trimmed.nc ../data/std30fldmean.nc ../data/anos_norm.nc

#mamba deactivate




