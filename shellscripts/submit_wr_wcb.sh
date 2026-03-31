#!/bin/bash
set -eo pipefail

START_YEAR=1979
END_YEAR=2024
MAX_CONCURRENT=8

mkdir -p logs

jid=$(sbatch --parsable --array=${START_YEAR}-${END_YEAR}%${MAX_CONCURRENT} wr_wcb_array.slurm)
echo "Submitted array job: ${jid}"

rid=$(sbatch --parsable --dependency=afterok:${jid} wr_wcb_reduce.slurm)
echo "Submitted reducer job: ${rid}"