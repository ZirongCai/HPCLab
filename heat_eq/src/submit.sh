#!/bin/bash
#SBATCH -o test.out
#SBATCH --time=00:15:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH -N 2
##SBATCH --cpus-per-task=1
##SBATCH --ntasks=2

module unload mpi.intel
module load openmpi

resolutions=(200 500 1200)

# Bash dosen't know 2d Arrays to just loop over two arrays
topox=( 1 3 7 4)
topoy=( 4 3 1 4)

for ((i=0;i<${#topox[@]};++i)); do
    for ((j=0;j<${#resolutions[@]};++j)); do
        ((num_procs=${topox[i]} * ${topoy[i]}))
        # printf "%d * %d = %d\n" ${topox[i]} ${topoy[i]} ${num_procs}
        printf "*************** %s x %s @ %s *********************\n" "${topox[i]}" "${topoy[i]}" "${resolutions[j]}"
        mpirun --oversubscribe -np ${num_procs} ./heat "input_${resolutions[j]}.dat" ${topox[i]} ${topoy[i]}
        printf "\n"
    done
    printf "\n\n\n"
done
