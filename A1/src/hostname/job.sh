#!/bin/bash
#SBATCH -J hello_host
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --export=NONE
#SBATCH --time=00:02:00
  
module load slurm_setup

echo "Running Program"
mpiexec -n $SLURM_NTASKS ./hello_host

