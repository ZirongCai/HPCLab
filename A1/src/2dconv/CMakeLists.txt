cmake_minimum_required(VERSION 2.8)

project(mm_examples)
set(CMAKE_CXX_STANDARD 11)

################################################################################
# build a helper(common) library
find_package(Qt5Widgets)

set(common_target "mm_simd")

add_executable(${common_target} common/helper.cpp 
                             common/compute_kernels_default.cpp
                             src/main.cpp)

target_compile_options(${common_target} PRIVATE "-O2")
target_include_directories(${common_target} PRIVATE src common)


add_subdirectory(simd)

target_link_libraries(${common_target} PRIVATE Qt5::Widgets -lboost_program_options conv_simd -lirc)

option(GRAPHICS "enable rendering of images" ON)
if (${GRAPHICS})
    target_compile_definitions(${common_target} PRIVATE GRAPHICS)
endif()
message(STATUS "graphics is: ${GRAPHICS}")


add_dependencies(${common_target} job_script)
add_custom_target(job_script
                 COMMAND cp ../scripts/job.sh .
                 COMMENT "Generating a SLURM job script")

