#!/bin/bash
#SBATCH -J cg_parallel
#SBATCH -D ./
#SBATCH --time=02:00:00
#SBATCH --mem=10000
#SBATCH --nodes=4
###SBATCH --ntasks-per-node=28
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
##SBATCH --clusters=cm2
##SBATCH --partition=cm2_std
##SBATCH --qos=cm2_std

module load slurm_setup

#ListOMP=(1 2 4 7 14) 
ListOMP=(1 7) 
nCores=$((${SLURM_JOB_NUM_NODES}*28))
RES="output_hybrid_${SLURM_JOB_NUM_NODES}node.txt"
TMP="out_hybrid${SLURM_JOB_NUM_NODES}.txt"

if [ ${SLURM_JOB_NUM_NODES} -eq 1 ]; then
  #ListTopoX=(7 7 7 2 2)
  ListTopoX=(7 2)
elif [ ${SLURM_JOB_NUM_NODES} -eq 2 ]; then
  #ListTopoX=(7 7 7 4 2)
  ListTopoX=(7 4)
elif [ ${SLURM_JOB_NUM_NODES} -eq 4 ]; then
  #ListTopoX=(7 7 7 4 4)
  ListTopoX=(7 4)
elif [ ${SLURM_JOB_NUM_NODES} -eq 8 ]; then
  #ListTopoX=(14 7 7 4 4)
  ListTopoX=(14 8)
fi

export KMP_AFFINITY=granularity=core,compact,1,0

printf "%s \t %s \t %s \t %s \t %s \t %s \t %s \n " "nCores" "nOMP" "TopoX" "TopoY" "nIterations" "Residum" "Time" > $RES

for i in "${!ListOMP[@]}"; do
  nOMP="${ListOMP[i]}"
  export OMP_NUM_THREADS=$nOMP
  nMPI=$(($nCores/$nOMP))
  nMPI_per_node=$(($nMPI/${SLURM_JOB_NUM_NODES}))
  nx="${ListTopoX[i]}"
  ny=$(($nMPI/$nx))
  export I_MPI_PIN_DOMAIN=omp I_MPI_PIN_ORDER=compact I_MPI_PIN_CELL=core I_MPI_PIN_CELL=core 
  srun -n ${nMPI} --mpi=pmi2 ./parallel 0.0001 5000 0.0005 $nx $ny > ${TMP}
  #----------------------------------
  time1=$(grep 'time' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  niter=$(grep 'Number of iterations' ${TMP} | grep -oE '[0-9]+ ')
  residum=$(grep 'Final norm of residuum' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  #---------------------------------
  printf "%s \t %s \t  %s \t %s \t %s \t %s \t %s\n" "$nCores" "$nOMP" "$nx" "$ny" "$niter" "$residum" "$time1"  >> $RES
  cat ${TMP}
  rm ${TMP}
done

