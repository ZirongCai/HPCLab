#!/bin/bash
#SBATCH -J Pi_StrongScale
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1-1
# 56 is the maximum reasonable value for CooLMUC-2
#SBATCH --mail-type=end
#SBATCH --mail-user=zirong.cai@tum.de
#SBATCH --export=NONE
#SBATCH --time=08:00:00
module load slurm_setup

export KMP_AFFINITY=balanced #scatter, unbalanced
printf "%s\t%s\t%s\n " "NumThread" "DataSize" "ExeTime" > $RES
#echo "nOMP     DataSize		ExeTime" > $RES

NumThreads="1 2 4 7 8 14 16 28 32 56"
ArrayLength="5000000 10000000 20000000 40000000"
for len in $ArrayLength; do
	RES="StrongScale_${len}.txt"
	for i in $NumThreads ; do
		export OMP_NUM_THREADS=$i
		./quicksort ${len} > out.txt
		printf "%s \t %s \t %s \n " "${i}" "${len}" $(tail -n1 out.txt | cut -d ":" -f 3 ) >> $RES
		rm out.txt	
	done
done