#!/bin/bash
export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK
case $OMPI_COMM_WORLD_LOCAL_RANK in
[0]) cpus=0;;
[1]) cpus=1;;
[2]) cpus=2;;
[3]) cpus=3;;
esac
numactl --physcpubind=$cpus $@
