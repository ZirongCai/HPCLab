#!/bin/bash
#SBATCH -J 2Dconv
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=$EMAIL
#SBATCH --export=ALL
#SBATCH --time=00:02:00


echo "Running default implementation"
./mm_simd --mode default -m 3 -k 3 

echo
echo
echo 

echo "Running SIMD implementation"
./mm_simd --mode intrinsics -m 3 -k 3
