#!/bin/bash
#SBATCH -J dgemm
#SBATCH -D ./
#SBATCH --time=04:00:00
#SBATCH --mem=8000
#SBATCH --nodes=1
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --mail-user=zirong.cai@tum.de


module load slurm_setup

KMP_AFFINITY=granularity=core,balanced,1,0
OMP_WAIT_POLICY=active
KMP_HOT_TEAMS=1
KMP_HOT_TEAMS_MAX_LEVELS=2
OMP_MAX_ACTIVE_LEVELS=2

RES="output_tasknum_testing.txt"
TMP="out28_2.txt"
printf "%s \t %s \t %s \t %s\n " "nTeams" "Microkernel" "GEBP" "GEMM" > $RES

for i in 1 2 4 7 14; do
  OMP_NUM_THREADS=28
  ./dgemm_task 2048 1 ${i} 100 > ${TMP}
  printf "%s \t %s \t %s \t %s\n" "${i}" $(cat ${TMP} | cut -d " " -f 4 | tr "\n" "\t") >> $RES
  rm ${TMP}
done
