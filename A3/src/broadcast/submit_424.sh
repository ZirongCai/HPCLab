#!/bin/bash
#SBATCH -J broadcast
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=zirong.cai@tum.de
#SBATCH --export=NONE
#SBATCH --time=08:00:00
  
module load slurm_setup


DataSize="2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152"
Num_Process="7 8 14 16 28 32"


for size in $DataSize; do

	RES="result_424_${size}.txt"
	echo "DataSize     nPRC    TrivialExeTime 	TreeExeTime 	MPI_BcastExeTime" > $RES
	for i in $Num_Process ; do
		mpiexec -n ${i} ./broadcast ${size} > out.txt
		printf "%s 		%s 		%s 		 %s 		 %s \n " "${size}" "${i}" $(tail -n3 out.txt | cut -d " " -f 6 ) | tr '\n' '\t'  >> $RES
		rm out.txt	
	done
done
