#!/bin/bash
#SBATCH
#SBATCH -J dgemm
#SBATCH -D ./
##SBATCH --get-user-env
##SBATCH -o result_compute.out_32x32
#SBATCH --time=02:00:00
#SBATCH --mem=8000
##SBATCH --clusters=serial
##SBATCH --partition=serial_std
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH -N 1
#SBATCH --cpus-per-task=1

module load slurm_setup

./dgemm
./dgemm
./dgemm
