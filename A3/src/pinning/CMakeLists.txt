cmake_minimum_required(VERSION 3.5)
project(pinning)


add_executable(${CMAKE_PROJECT_NAME} main.cpp)
target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE -fopenmp)
target_link_libraries(${CMAKE_PROJECT_NAME} -fopenmp)
