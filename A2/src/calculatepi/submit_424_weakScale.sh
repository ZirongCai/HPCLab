#!/bin/bash
#SBATCH -o test_421.out
#SBATCH --time=00:05:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1
#SBATCH --cpus-per-task=48


affinities="compact scatter balanced"
NumThreads="1 2 4 8 14 16 28 32 56"

for aff in $affinities; do
	export KMP_AFFINITY=$aff
	RES="test${aff}.txt"
	echo "NumThread		Precision		SerialExeTime		CriticalExeTime		ReductionExeTime" > $RES

	for i in $NumThreads ; do
		precision=$(expr $i \* 1000)
		export OMP_NUM_THREADS=$i
		./calculate_pi $precision > out.txt
		printf "%s \t  %s \t  %s  \t  %s \t  %s " "${i}" "${precision}" $(tail -n3 out.txt | cut -d " " -f 6 ) | tr '\n' '\t' >> $RES
		printf "\n" >> $RES
		rm out.txt	
	done
done
