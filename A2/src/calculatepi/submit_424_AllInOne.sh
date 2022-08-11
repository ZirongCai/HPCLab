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



NumThreads="1 2 4 7 8 14 16 28 56"
affinities="compact scatter balanced"
PRECISION="1000 2000 4000 7000 8000 14000 16000 28000 56000"

for aff in $affinities; do
	export KMP_AFFINITY=$aff
	for precision in $PRECISION ; do
		RES="scaling${aff}_${precision}.txt"
		echo "NumThread		Precision		SerialExeTime		CriticalExeTime		ReductionExeTime" > $RES
		for i in $NumThreads ; do
			export OMP_NUM_THREADS=$i
			./calculate_pi ${precision} > out.txt
			printf "%s \t  %s \t  %s  \t  %s \t  %s " "${i}" "${precision}" $(tail -n3 out.txt | cut -d " " -f 6 ) | tr '\n' '\t' >> $RES
			printf "\n" >> $RES
			rm out.txt	
		done
	done
done
