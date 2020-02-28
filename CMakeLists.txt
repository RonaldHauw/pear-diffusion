cmake_minimum_required(VERSION 3.15)
project(pear_diffusion)

set(CMAKE_CXX_STANDARD 14)

add_executable(pear_diffusion_2 main.cpp vector.hpp)

add_compile_options(-Wall)