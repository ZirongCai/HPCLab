#!/bin/bash

NumMPI="1 2 4"
GridSize="512 1024 2048 4096 8192"
#TODO
export NSTEPS=1
# end                     
#end

spack load openmpi/cc
spack load cuda/zi

#nvidia-smi
#mpirun --version
nvidia-smi -pm 1

for size in $GridSize; do
    export GRID_SIZE_X=${size}
    export GRID_SIZE_Y=${size}
	for mpi in $NumMPI ; do
        export nMPI=${mpi}
        APP="./build/swe-perf -x ${GRID_SIZE_X} -y $GRID_SIZE_Y -c $NSTEPS -o ."
		mpirun -np ${nMPI} ./gpu_bind.sh nsys profile -o ./nsys_profile/cuda-aware/profile_cuda_aware_${size}_${nMPI} --trace=mpi,cuda --mpi-impl openmpi $APP
	done
done

echo "DONE\n"
