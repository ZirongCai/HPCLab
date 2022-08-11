cmake_minimum_required(VERSION 2.8)


#include(CMakeForceCompiler)
#CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
#CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")

add_library(conv_simd convolution.cpp)

target_include_directories(conv_simd PUBLIC)
target_compile_options(conv_simd PRIVATE "-O2" "-march=native")
target_include_directories(conv_simd PRIVATE ../src ../common)
