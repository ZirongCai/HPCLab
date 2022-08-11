# How to compile
To compile our best CUDA-Aware MPI version, run as follow:  
`$ cmake -DENABLE_CUDA=ON -DENABLE_CUDA_AWARE_MPI=ON ..`  
`$ make -j`  

The CMakeFiles.txt provide the following options:  
`-DENABLE_TEST` (=OFF by default): Compile the code which has TEST feature (printing out WaterHeight).    
`-DENABLE_OUTPUT` (=ON by default): output vts files, one might turn it off for performance tests.    
`-DENABLE_NON_BLOCKING_COMM` (=ON by default): non-blocking communication for CUDA-Aware MPI implementation.  
`-DENABLE_CUDA_MEM_ADVISE` (=ON by default): enable MemAdvise for CUDA-Aware MPI.  
`-DENABLE_PACKING` (=ON by default): enable packing data for MPI communications in CUDA-Aware MPI.  
`-DENABLE_CUDA_INIT_OPT` (=OFF by default): "First touch" implementation for CUDA-Aware MPI.  

Important note: one always needs to have -DENABLE_CUDA=ON when compile the CUDA_AWARE_MPI.  

# How to run

### Normal simulation
One can uses the gpu_bind.sh script to execute multiple GPUs, as follows:  
`$ mpirun -np 2 ./gpu_bind.sh ./build/swe-mpi -x 4096 -y 4096 -c 4 -o outputnames_`

### Error test
First, set the grid size parameter in test/run_test.sh, the default test is with NX 32 NY 32 COUNT 5.
Then execute the test script with:
`$ bash run_test.sh`

### Profiling
`$ bash nsys_profile/profile_opt.sh`

