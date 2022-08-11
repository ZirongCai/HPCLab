#!/bin/bash
#SBATCH -J cg_parallel
#SBATCH -D ./
#SBATCH --time=04:00:00
#SBATCH --mem=8000
#SBATCH --nodes=1
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
##SBATCH --hint=nomultithread
#SBATCH --ntasks-per-node=28

module load slurm_setup

RES="output_gridsizes.txt"
TMP="out_gridsizes.txt"

nMPI=28

printf "%s \t %s \t %s \t %s \n " "Gridsize" "nIterations" "Residum" "Time" > $RES
grid=0.1
for i in {1..10}; do
  ###################
  grid=$(perl -E "say $grid/2")
  I_MPI_PIN_CELL=core mpiexec -n $SLURM_NTASKS ./parallel $grid 5000 0.0005 4 7 > ${TMP}
  #-----------------
  time1=$(grep 'time' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  niter=$(grep 'Number of iterations' ${TMP} | grep -oE '[0-9]+ ')
  residum=$(grep 'Final norm of residuum' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  printf "%s \t %s \t %s \t %s\n" "$grid" "$niter" "$residum" "$time1"  >> $RES
  #-----------------
  cat ${TMP}
  rm ${TMP}
done

