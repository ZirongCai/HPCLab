#!/bin/bash

#module unload intel-mpi
#module unload intel/19.0
 
#module load gcc
#module load mpi.ompi/2.1/gcc

module load boost
module load qt

export CPATH=$BOOST_INCDIR:$CPATH
export C_INCLUDE_PATH=$BOOST_INCDIR:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$BOOST_INCDIR:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$BOOST_LIBDIR:$LIBRARY_PATH

