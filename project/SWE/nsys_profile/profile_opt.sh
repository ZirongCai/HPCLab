#!/bin/bash

#NumMPI="1 2 4"
#GridSize="512 1024 2048 4096 8192"
NumMPI="4"
GridSize="512"
#TODO
export NSTEPS=1
save_to="output_opt"
# end                     

spack load openmpi/cc
spack load cuda/zi

nvidia-smi -pm 1
mkdir -p $save_to

for size in $GridSize; do
    export GRID_SIZE_X=${size}
    export GRID_SIZE_Y=${size}
	for mpi in $NumMPI ; do
        export nMPI=${mpi}
        APP="../build/swe-mpi -x ${GRID_SIZE_X} -y $GRID_SIZE_Y -c $NSTEPS -o ."
		mpirun -np ${nMPI} ../gpu_bind.sh nsys profile -o profile_cuda_aware_${size}_${nMPI} --trace=mpi,cuda --mpi-impl openmpi $APP
        nsys stats profile_cuda_aware_${size}_${nMPI}.qdrep > $save_to/statsout_cuda_aware_${size}_${nMPI}.txt
	done
done

rm ._*.nc *.qdrep

echo "DONE\n"
