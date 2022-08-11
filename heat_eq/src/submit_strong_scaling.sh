#!/bin/bash
#SBATCH -o test.out
#SBATCH --time=00:05:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=48  

module unload mpi.intel
module load openmpi

topox=(1 2 3 4 6 8 12 16 24 48)
topoy=(48 24 16 12 8 6 4 3 2 1)

nNode=$SLURM_JOB_NUM_NODES
  
re_time="../result_logs/re71_strongscale_${nNode}nodes_time.txt"
re_perf="../result_logs/re71_strongscale_${nNode}nodes_perf.txt"
re_all="../result_logs/re71_strongscale_${nNode}nodes_alldata.txt"

rm -f $re_time $re_perf $re_all

max_cores=$(($nNode * 48))
 
#n=$(($max_cores))
#for i in $(seq 1 $n)
#do
# [ $(expr $n / $i \* $i) == $n ] && echo $i $(($n / $i))
#done

for i in "${!topox[@]}"; do 
    x=$((${topox[i]})) 
    y=$((${topoy[i]}*$nNode)) 
    nproc=$(($x*$y))
    
    echo $nproc
     
    # Exit if we go over max cores 
    if [ $nproc -gt $max_cores ]; then
        continue
    fi

    echo "mpirun -n $nproc heat small.dat $x $y"
    mpirun -np $nproc ./heat test.dat $x $y 2>&1 | tee out.txt

    # Extract data for each each, put time to re_time, GFLOPS to re_perf, and save raw data to re_all	
    echo "---------------Topology: ${x} x ${y} ------------------" | tee -a $re_time $re_perf $re_all > /dev/null
    printf "[%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s] \n" $(tail -n11 out.txt | awk '{print $4}'| tr '\n' ' ') >> $re_time
    printf "[%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s] \n" $(tail -n11 out.txt | awk '{print $8}'| tr '\n' ' ') >> $re_perf
    cat out.txt >> $re_all
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $re_all
    echo " " | tee -a $re_time $re_perf $re_all > /dev/null
    rm out.txt 
done

## This script is for strong-scalling test for assisgment 7.1, change input file to test.dat, and proper time (max 15 mins on testing partition)
