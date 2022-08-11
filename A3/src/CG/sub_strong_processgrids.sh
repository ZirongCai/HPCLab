#!/bin/bash
#SBATCH -J cg_parallel
#SBATCH -D ./
#SBATCH --time=04:00:00
#SBATCH --mem=10000
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=28
##SBATCH --clusters=cm2_tiny
##SBATCH --partition=cm2_tiny
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --qos=cm2_std

module load slurm_setup

RES="output_strong_8nodes.txt"
TMP="out_strong8.txt"

printf "%s \t %s \t %s \t %s \t %s \t %s \n " "nMPI" "TopoX" "TopoY" "nIterations" "Residum" "Time" > $RES
nMPI=$SLURM_NTASKS

for i in 1 2 4 7 14 28 56; do
  j=$(($nMPI/$i))
  I_MPI_PIN_CELL=core mpiexec -n $SLURM_NTASKS ./parallel 0.00025 5000 0.0005 ${i} ${j} > ${TMP}
  #----------------------------------
  time1=$(grep 'time' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  niter=$(grep 'Number of iterations' ${TMP} | grep -oE '[0-9]+ ')
  residum=$(grep 'Final norm of residuum' ${TMP} | grep -oE "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.")
  #---------------------------------
  printf "%s \t %s \t %s \t %s \t %s \t %s\n" "$nMPI" "$i" "$j" "$niter" "$residum" "$time1"  >> $RES
  cat ${TMP}
  rm ${TMP}
done

