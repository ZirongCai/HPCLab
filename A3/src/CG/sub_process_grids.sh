#!/bin/bash
#SBATCH -J cg_parallel
#SBATCH -D ./
#SBATCH --time=04:00:00
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
##SBATCH --hint=nomultithread

module load slurm_setup

RES="output_processgrids.txt"
TMP="out_process.txt"

nMPI=28

printf "%s \t %s \t %s \t %s \t %s \t %s \n " "nMPI" "TopoX" "TopoY" "nIterations" "Residum" "Time" > $RES

for i in 1 2 4 7 14 28; do
  j=$(($nMPI/$i))
  I_MPI_PIN_CELL=core mpiexec -n $SLURM_NTASKS ./parallel 0.0005 5000 0.0005 ${i} ${j} > ${TMP}
  #----------------------------------
  time1=$(grep 'time' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  niter=$(grep 'Number of iterations' ${TMP} | grep -oE '[0-9]+ ')
  residum=$(grep 'Final norm of residuum' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  #---------------------------------
  printf "%s \t %s \t %s \t %s \t %s \t %s\n" "$nMPI" "$i" "$j" "$niter" "$residum" "$time1"  >> $RES
  cat ${TMP}
  rm ${TMP}
done

