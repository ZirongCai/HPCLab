#!/bin/bash
#SBATCH -o test_421.out
#SBATCH --time=00:05:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1
#SBATCH --cpus-per-task=48


NumThreads="1 2 4 8 12 16 24 32 48"
affinities="compact scatter balanced"

for aff in $affinities; do
	export KMP_AFFINITY=$aff
	RES="result_424_res6500_${aff}.txt"
	echo "nOMP     ExeTime         Residual        MFLOPS     Flop_instruction(M)" > $RES

	for i in $NumThreads ; do
		export OMP_NUM_THREADS=$i
		./heat test.dat > out.txt
		printf "%s \t %s \t %s \t %s \t %s \n " "${i}" $(tail -n5 out.txt | cut -d ":" -f 2 | tr '\n' ' ') >> $RES
		rm out.txt	
	done
done
