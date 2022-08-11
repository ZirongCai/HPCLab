#!/bin/bash
#SBATCH -J cg_scorep
#SBATCH -D ./
#SBATCH --time=01:00:00
#SBATCH --mem=10000
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=28
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny

module load slurm_setup
module load cube scorep scalasca

export OMP_NUM_THREADS=7
export KMP_AFFINITY=granularity=core,compact,1,0
export I_MPI_PIN_DOMAIN=omp I_MPI_PIN_ORDER=compact I_MPI_PIN_CELL=core I_MPI_PIN_CELL=core
scalasca -instrument mpiexec -np 4 ./parallel 0.00025 5000 0.0005 2 2
