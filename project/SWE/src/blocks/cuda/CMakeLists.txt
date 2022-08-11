#Todo: add other solvers

if (CUDA_AWARE_MPI)
  set (SELECTED_SWE_BLOCK_CUDA SWE_BlockCudaAwareMPI.cu)
  set (SELECTED_SWE_BLOCK_CUDA_HEADER SWE_BlockCudaAwareMPI.hh)

  set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA SWE_WavePropagationBlockCudaAwareMPI.cu)
  set(SELECTED_SWE_BLOCK_HEADER_CUDA SWE_WavePropagationBlockCudaAwareMPI.hh)

  set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.cu)
  set(SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.hh)

else()
  set (SELECTED_SWE_BLOCK_CUDA SWE_BlockCUDA.cu)
  set (SELECTED_SWE_BLOCK_CUDA_HEADER SWE_BlockCUDA.hh)

  if (SOLVER STREQUAL "fwave")
    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA SWE_WavePropagationBlockCuda.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA SWE_WavePropagationBlockCuda.hh)

    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.hh)
  elseif(SOLVER STREQUAL "rusanov")
    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA SWE_RusanovBlockCUDA.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA SWE_RusanovBlockCUDA.hh)

    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS SWE_RusanovBlockCUDA_kernels.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS SWE_RusanovBlockCUDA_kernels.hh)
  elseif(SOLVER STREQUAL "augRie")
    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA SWE_WavePropagationBlockCuda.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA SWE_WavePropagationBlockCuda.hh)

    set(SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.cu)
    set(SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS SWE_WavePropagationBlockCuda_kernels.hh)
  else()
    message (FATAL_ERROR "${SOLVER} is not supported with cuda")
  endif()
endif()

target_sources(swe 
    PRIVATE
    ${SELECTED_SWE_BLOCK_CUDA}
    SWE_BlockCUDA_kernels.cu
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA}
    ${SELECTED_SWE_BLOCK_COMPILATION_UNIT_CUDA_KERNELS}
    PUBLIC
    ${SELECTED_SWE_BLOCK_CUDA_HEADER}
    SWE_BlockCUDA_kernels.hh
    ${SELECTED_SWE_BLOCK_HEADER_CUDA}
    ${SELECTED_SWE_BLOCK_HEADER_CUDA_KERNELS}
)
