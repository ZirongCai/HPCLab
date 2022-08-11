#!/bin/bash


NumMPI="1 2 4"
GridSize="64 128 256 512 1024 2048 4096 8192"
#TODO
export NSTEPS=1
# end                     
#end

spack load openmpi/cc
spack load cuda/zi

#nvidia-smi
#mpirun --version
nvidia-smi -pm 1
mkdir -p outputs

for size in $GridSize; do
    export GRID_SIZE_X=${size}
    export GRID_SIZE_Y=${size}
    RES="StrongScale_CUDA_${size}.txt"
	for mpi in $NumMPI ; do
        export nMPI=${mpi}
		mpirun -np ${nMPI} -bind-to none ./profile/execute_swe_${mpi}gpu.sh >> ./profile/outputs/$RES
	done
done

echo "DONE\n"
