cmake_minimum_required(VERSION 3.15)
project(pear_diffusion)

set(CMAKE_CXX_STANDARD 14)

SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -Wall -fstack-protector")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(pear_diffusion software/main.cpp)

add_compile_options(-Wall)
add_compile_options(-O3)
add_compile_options(-fstack-protector)
