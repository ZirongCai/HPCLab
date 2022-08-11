#!/bin/bash
#SBATCH -J dgemm
#SBATCH -D ./
#SBATCH --time=04:00:00
#SBATCH --mem=8000
#SBATCH --nodes=1
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny

module load slurm_setup

KMP_AFFINITY=granularity=core,balanced,1,0
OMP_WAIT_POLICY=active
KMP_HOT_TEAMS=1
KMP_HOT_TEAMS_MAX_LEVELS=2
OMP_MAX_ACTIVE_LEVELS=2

RES="output_weak.txt"
TMP="outw${size}.txt"
printf "%s \t %s \t %s\n " "nThreads" "Size" "Time" > $RES

nThreads=("1" "2" "4" "7" "14" "28")
sizes=("512" "640" "800" "960" "1216" "1536")

for i in {0..3}; do
  OMP_NUM_THREADS=${nThreads[${i}]}
  ./dgemm ${sizes[${i}]} 1 ${nThreads[${i}]} 100 > ${TMP}
  time=$(cat ${TMP} | tail -n1 | cut -d " " -f 2)
  printf "%s \t %s \t %s\n" "${nThreads[${i}]}" "${sizes[${i}]}" "${time}" >> $RES
  rm ${TMP}
done
for i in {4..5}; do
  OMP_NUM_THREADS=${nThreads[${i}]}
  ./dgemm ${sizes[${i}]} 1 2 100 > ${TMP}
  time=$(cat ${TMP} | tail -n1 | cut -d " " -f 2)
  printf "%s \t %s \t %s\n" "${nThreads[${i}]}" "${sizes[${i}]}" "${time}" >> $RES
  rm ${TMP}
done
