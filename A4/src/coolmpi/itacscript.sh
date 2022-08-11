#/bin/bash
module unload intel-mpi intel-mkl intel
module load intel-parallel-studio/2018
#everything is now in $INTEL_PARALLEL_STUDIO_BASE
. $INTEL_PARALLEL_STUDIO_BASE/bin/compilervars.sh intel64
. $INTEL_PARALLEL_STUDIO_BASE/vtune_amplifier/amplxe-vars.sh
. $INTEL_PARALLEL_STUDIO_BASE/itac_latest/bin/itacvars.sh

#mpirun -n X -trace ./cool_mpi_app_exe a b c1 c2 c3 c4 c5 cn

