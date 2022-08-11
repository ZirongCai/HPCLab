#!/bin/bash

# if you can not compile the code on supermuc (it raises error: can not find mpi.h), source this file and try again!

module unload mpi.intel
module load openmpi

export "CPATH=/dss/dsshome1/lrz/sys/spack/release/19.2/opt/x86_avx512/openmpi/4.0.2-intel-25wtmyq/include:$CPATH"
export CC=icc 
