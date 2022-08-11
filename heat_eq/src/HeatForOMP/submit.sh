#!/bin/bash
#SBATCH -o result_sequential.txt
#SBATCH --time=00:15:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 1

module unload mpi.intel
module load openmpi

./heat test.dat
