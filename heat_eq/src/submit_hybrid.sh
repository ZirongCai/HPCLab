#!/bin/bash
#SBATCH -o test.out
#SBATCH --time=00:15:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1

module unload mpi.intel
module load openmpi

nNode=${SLURM_JOB_NUM_NODES}

Threads="1 2 4 8 12 24 48"

#re_time="../result_logs/hybrid_${nnode}node_time.txt"
#re_perf="../result_logs/hybrid_${nnode}node_perf.txt"
#re_all="../result_logs/hybrid_${nnode}node_alldata.txt"

re_time="save_results/hybrid_${nNode}node_time.txt"
re_perf="save_results/hybrid_${nNode}node_perf.txt"
re_all="save_results/hybrid_${nNode}node_alldata.txt"
rm -f $re_time $re_perf $re_all 

for nOMP in $Threads; do
	nMPI=$(($nNode*48/$nOMP)) #with Threads list above, 1-4 nodes, possible nMPI is 8 16 24 48 96
     	export OMP_NUM_THREADS=$nOMP
    
	if [[ $nMPI -eq 1 ]]
	then
	  topox="1"
	elif [[ $nMPI -eq 2 ]]
	then
	  topox="1 2"
	elif [[ $nMPI -eq 3 ]]
	then
	  topox="1 2 3"
	elif [[ $nMPI -eq 4 ]]
	then
	  topox="1 2 4"
	elif [[ $nMPI -eq 6 ]]
	then
	  topox="1 2 3 6"
	elif [[ $nMPI -eq 8 ]]
	then
	  topox="1 2 4 8"
	elif [[ $nMPI -eq 12 ]]
	then
	  topox="1 2 3 6 12"
	else	
	  topox="1 2 4 8 12 24"
	fi

	echo "*** ${nMPI} MPI and ${nOMP} OMP ******************************" 2>&1 | tee -a $re_time $re_perf $re_all

	for x in $topox; do
		y=$(($nMPI/$x))
         
		echo "---------------Topology: ${x} x ${y} ------------------" 2>&1 | tee -a $re_time $re_perf $re_all
		
		srun -n $nMPI heat test.dat $x $y > out.txt
 
		# Extract data for each each, put time to re_time, GFLOPS to re_perf, and save raw data to re_all	
		printf "[%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s]\n" $(tail -n11 out.txt | awk '{print $4}'| tr '\n' ' ') >> $re_time
		printf "[%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s]\n" $(tail -n11 out.txt | awk '{print $8}'| tr '\n' ' ') >> $re_perf
		cat out.txt >> $re_all
		echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $re_all
		echo " " | tee -a $re_time $re_perf $re_all > /dev/null
		#rm out.txt	
	done
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" >> $re_all
    echo "" 
    echo "" 
    echo "" 
    echo "" 
done

## This script is for hybrid test in ass 7.2
