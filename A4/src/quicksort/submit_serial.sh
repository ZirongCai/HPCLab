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

#export KMP_AFFINITY=compact #scatter, unbalanced
printf "%s\t%s\t%s\n " "NumThread" "DataSize" "ExeTime" > $RES
#echo "nOMP     DataSize		ExeTime" > $RES
RES="StrongScale_Serial.txt"

ArrayLength="5000000 10000000 20000000 40000000"
for len in $ArrayLength; do
	./quicksort ${len} > out.txt
	printf "%s \t %s \t %s \n " "${i}" "${len}" $(tail -n1 out.txt | cut -d ":" -f 3 ) >> $RES
	rm out.txt	
done