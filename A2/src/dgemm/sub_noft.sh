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

size=2048
RES="outputnew_${size}s_noft.txt"
TMP="out${size}_no.txt"
printf "%s \t %s \t %s \t %s\n " "nThreads" "Microkernel" "GEBP" "GEMM" > $RES

for i in 1 2 4 7; do
  OMP_NUM_THREADS=${i}
  ./dgemm_no ${size} 1 ${i} 500 > ${TMP}
  printf "%s \t %s \t %s \t %s\n" "${i}" $(cat ${TMP} | cut -d " " -f 4 | tr "\n" "\t") >> $RES
  echo "******* nThreads ${i} ********"
  cat ${TMP}
  rm ${TMP}	
done
for i in 14 28 56; do
  OMP_NUM_THREADS=${i}
  ./dgemm_no ${size} 1 2 500 > ${TMP}
  printf "%s \t %s \t %s \t %s\n" "${i}" $(cat ${TMP} | cut -d " " -f 4 | tr "\n" "\t") >> $RES
  echo "******* nThreads ${i} ********"
  cat ${TMP}
  rm ${TMP}	
done
