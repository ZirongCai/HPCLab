#!/bin/bash
#SBATCH -o test_421.out
#SBATCH --time=00:05:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1
#SBATCH --cpus-per-task=48

#export KMP_AFFINITY=compact #scatter, unbalanced
RES="StrongScale.txt"
printf "%s \t %s \t %s \t %s \t %s \n " "NumThread" "Precision" "SerialExeTime" "CriticalExeTime" "ReductionExeTime" > $RES
#echo "nOMP     DataSize		ExeTime" > $RES

NumThreads="1 2 4 7 14 28 56"
Precision="1000 2000 4000 8000 16000 32000"
for i in $NumThreads ; do
	export OMP_NUM_THREADS=$i
	for j in $Precision ; do
		./calculate_pi $j > out.txt
		printf "%s \t  %s \t  %s  \t  %s \t  %s " "${i}" "${j}" $(tail -n3 out.txt | cut -d " " -f 6 ) | tr '\n' '\t' >> $RES
		printf "\n" >> $RES
		rm out.txt	
	done
done
