# HPC Lab

- The SWE project is in main branch inside project/SWE directory.

- To run the code
```
cd SWE
mkdir build
cd build
cmake -DENABLE_CUDA=ON -DENABLE_CUDA_AWARE_MPI=ON -D(Option listed in CMakelist) ..
then you should have the executable binary
```

- To test the code
```
cd SWE/test
bash ./run_test.sh
```
