# Parameter
NX=32
NY=32
STEPS=5

## Build the binaries
spack load cuda/zif openmpi/ccn netcdf-c@4.6.1
mkdir -p build
cd build

if [ ! -f swe-no-cuda ]; then
cmake -DENABLE_CUDA=OFF -DENABLE_CUDA_AWARE_MPI=OFF -DENABLE_ERROR_TEST=ON -DENABLE_OUTPUT=OFF ../../
make -j
cp swe-mpi swe-no-cuda
make clean
fi

if [ ! -f swe-cuda ]; then
cmake -DENABLE_CUDA=OFF -DENABLE_CUDA_AWARE_MPI=OFF -DENABLE_ERROR_TEST=ON -DENABLE_OUTPUT=OFF ../../
make -j
cp swe-mpi swe-cuda
make clean
fi

if [ ! -f swe-cuda-aware ]; then
cmake -DENABLE_CUDA=OFF -DENABLE_CUDA_AWARE_MPI=ON -DENABLE_ERROR_TEST=ON -DENABLE_OUTPUT=OFF ../../
make -j
cp swe-mpi swe-cuda-aware
make clean
fi

## Run all tests
mpirun -np 1 ./swe-no-cuda -x ${NX} -y ${NY} -c ${STEPS} -o nocuda_1mpi_ > stdout_1mpi.txt 2&>1
od -v -An -f WaterHeight_1mpi.txt | tr '\n' ' ' > out_1mpi.txt
#
mpirun -np 4 ./swe-no-cuda -x ${NX} -y ${NY} -c ${STEPS} -o nocuda_4mpi_ > stdout_4mpi.txt 2&>1
od -v -An -f WaterHeight_4mpi.txt | tr '\n' ' ' > out_4mpi.txt
#
mpirun -np 1 ./swe-cuda -x ${NX} -y ${NY} -c ${STEPS} -o cuda_1mpi_ > stdout_1cuda.txt 2&>1
od -v -An -f WaterHeight_1mpi.txt | tr '\n' ' ' > out_1cuda.txt
#
mpirun -np 4 ../gpu_bind.sh ./swe-cuda -x ${NX} -y ${NY} -c ${STEPS} -o cuda_4mpi_ > stdout_4cuda.txt 2&>1
od -v -An -f WaterHeight_4mpi.txt | tr '\n' ' ' > out_4cuda.txt
#
mpirun -np 1 ./swe-cuda-aware -x ${NX} -y ${NY} -c ${STEPS} -o cudaaware_1mpi_ > stdout_1cudaAware.txt 2&>1
od -v -An -f WaterHeight_1mpi.txt | tr '\n' ' ' > out_1cudaAware.txt
#
mpirun -np 4 ../gpu_bind.sh ./swe-cuda-aware -x ${NX} -y ${NY} -c ${STEPS} -o cudaaware_4mpi_ > stdout_4cudaAware.txt 2&>1
od -v -An -f WaterHeight_4mpi.txt | tr '\n' ' ' > out_4cudaAware.txt

## Error checking
echo "* TEST: 1 mpi vs. 1 cuda *"
python3 ../error_check.py out_1mpi.txt out_1cuda.txt

echo "* TEST: 1 mpi vs. 1 cuda-aware *"
python3 ../error_check.py out_1mpi.txt out_1cudaAware.txt

echo "* TEST: 1 cuda vs 1 cuda-aware *"
python3 ../error_check.py out_1cuda.txt out_1cudaAware.txt

echo "* TEST: 1 mpi vs. 4 mpi *"
python3 ../error_check.py out_1mpi.txt out_4mpi.txt

echo "* TEST: 1 cuda vs. 4 cuda *"
python3 ../error_check.py out_1cuda.txt out_4cuda.txt

echo "* TEST: 1 cuda-aware vs. 4 cuda-aware *"
python3 ../error_check.py out_1cudaAware.txt out_4cudaAware.txt

echo "* TEST: 1 mpi vs. 4 cuda *"
python3 ../error_check.py out_1mpi.txt out_4cuda.txt

echo "* TEST: 1 mpi vs. 4 cuda-aware *"
python3 ../error_check.py out_1mpi.txt out_4cudaAware.txt

echo "* TEST: 1 cuda vs. 4 mpi *"
python3 ../error_check.py out_1cuda.txt out_4mpi.txt

echo "* TEST: 1 cuda vs. 4 cuda-aware *"
python3 ../error_check.py out_1cuda.txt out_4cudaAware.txt

echo "* TEST: 1 cuda-aware vs. 4 cuda *"
python3 ../error_check.py out_1cudaAware.txt out_4cuda.txt

echo "* TEST: 4 mpi vs. 4 cuda *"
python3 ../error_check.py out_4mpi.txt out_4cuda.txt

echo "* TEST: 4 cuda vs 4 cuda-aware *"
python3 ../error_check.py out_4cuda.txt out_4cudaAware.txt

echo "* TEST: 4 mpi vs. 4 cuda-aware *"
python3 ../error_check.py out_4mpi.txt out_4cudaAware.txt

## Cleanup
rm *.txt
rm *.nc
cd ..
