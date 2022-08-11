#!/bin/bash
#SBATCH -J cg_parallel
#SBATCH -D ./
#SBATCH --time=04:00:00
#SBATCH --mem=10000
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=28
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
##SBATCH --hint=nomultithread

module load slurm_setup

RES="output_strong_new.txt"
TMP="out_strong.txt"

nMPI=(1 2 4 8 14 28)
TopoX=(1 1 2 2 2 4)

printf "%s \t %s \t %s \t %s \t %s \t %s \n " "nMPI" "TopoX" "TopoY" "nIterations" "Residum" "Time" > $RES

for i in "${!nMPI[@]}"; do
  nx="${TopoX[i]}"
  np="${nMPI[i]}"
  ny=$(($np/$nx))
  I_MPI_PIN_CELL=core mpiexec -n $np ./parallel 0.00025 5000 0.0005 $nx $ny > ${TMP}
  #----------------------------------
  time1=$(grep 'time' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  niter=$(grep 'Number of iterations' ${TMP} | grep -oE '[0-9]+ ')
  residum=$(grep 'Final norm of residuum' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  #---------------------------------
  printf "%s \t %s \t %s \t %s \t %s \t %s\n" "$np" "$nx" "$ny" "$niter" "$residum" "$time1"  >> $RES
  cat ${TMP}
  rm ${TMP}
done

