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

# <<'COMMENT'

ERA5_DIR=/daten/reana/arch/reanalysis/reanalysis/DKRZ/IFS/ERA5/day/atmos/zg/r1i1p1-050000Pa/

echo "Cut spatially and merge temporally"

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

echo "Calculating moving window std"

cdo runstd,30 temporally_merged.nc  std30.nc
cdo settime,12:00 std30.nc std30_12h.nc
cdo fldmean std30_12h.nc std30fldmean.nc

echo "Removing seasonal cycle by normalization with std"

cdo seldate,1950-01-16,2023-12-17 temporally_merged.nc trimmed.nc

cdo div trimmed.nc -enlarge,trimmed.nc std30fldmean.nc ../data/noSeasCycle.nc

rm *.nc

cdo -ydaymean ../data/noSeasCycle.nc ../data/ydailymean.nc

# COMMENT

echo "Calulating anomalies"

for Y in {1950..2023}
do
    echo "Year=$Y"
    cdo -sub -selyear,${Y} ../data/noSeasCycle.nc ../data/ydailymean.nc ../data/${Y}_ano.nc
done

cdo -O -mergetime ../data/*_ano.nc ../data/anos.nc

rm -f ../data/*_ano.nc

cdo -P 48 -detrend ../data/anos.nc ../data/anosdetrend.nc

# COMMENT



