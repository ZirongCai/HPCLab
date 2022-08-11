#!/bin/bash
#SBATCH
#SBATCH -J pinning
#SBATCH -D ./
#SBATCH --time=02:00:00
#SBATCH --mem=8000
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH -N 1


module load slurm_setup
echo "*** Task 1 *********************************************************"
OMP_NUM_THREADS=4 I_MPI_PIN_DOMAIN=omp I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core mpirun -np 2 ./pinning

echo "*** Task 2 *********************************************************"
OMP_NUM_THREADS=2 KMP_AFFINITY="explicit,proclist=[0,6,7,13,14,20,21,27]" I_MPI_PIN_ORDER=spread mpiexec -n 4 ./pinning 2> /dev/null

echo "*** Task 3 *********************************************************"
srun -n 2 --mpi=pmi2 --cpu-bind=map\_cpu:0,14  --hint=nomultithread ./pinning
echo ""
srun -n 4 --mpi=pmi2 --cpu-bind=map\_cpu:0,7,14,21  --hint=nomultithread ./pinning
