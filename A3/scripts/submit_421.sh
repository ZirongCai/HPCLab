#!/bin/bash
#SBATCH -o test_421.out
#SBATCH --time=00:05:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1
#SBATCH --cpus-per-task=48

#export KMP_AFFINITY=compact #scatter, unbalanced
RES="result_421_res6500.txt"
echo "nOMP     ExeTime         Residual        MFLOPS     Flop_instruction(M)" > $RES

NumThreads="1 2 4 8 12 16 24 32 48"

for i in $NumThreads ; do
	export OMP_NUM_THREADS=$i
	./heat test.dat > out.txt
	printf "%s \t %s \t %s \t %s \t %s \n " "${i}" $(tail -n5 out.txt | cut -d ":" -f 2 | tr '\n' ' ') >> $RES
	rm out.txt	
done
