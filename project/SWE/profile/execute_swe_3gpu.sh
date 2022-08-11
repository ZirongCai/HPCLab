#!/bin/bash

APP="./build/swe-perf -x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -c ${NSTEPS} -o ."

CPU_CORES_PER_RANK=1
lrank=$OMPI_COMM_WORLD_LOCAL_RANK

case ${lrank} in
0)
#ldd $APP

  #set GPU and CPU affinity of local rank 
  export CUDA_VISIBLE_DEVICES=0; numactl --cpunodebind=0  $APP
  ;;
1)

  #set GPU and CPU affinity of local rank 
  export CUDA_VISIBLE_DEVICES=1; numactl --cpunodebind=0  $APP

  ;;
2)
  #set GPU and CPU affinity of local rank 
  export CUDA_VISIBLE_DEVICES=2; numactl --cpunodebind=0  $APP

  ;;
esac
