cmake_minimum_required(VERSION 3.15)
project(pear_diffusion)

set(CMAKE_CXX_STANDARD 14)
SET(GCC_COVERAGE_COMPILE_FLAGS "-O1 -Wall -fsanitize=address")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(pear_diffusion_2 main.cpp)

add_compile_options(-Wall -fsanitize=address)