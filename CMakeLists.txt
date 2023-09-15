cmake_minimum_required(VERSION 3.1...3.23)
project(ga_spgmt_uu)

set(CMAKE_CXX_STANDARD 17)

# dependencies
add_subdirectory(dep/Catch2)
add_subdirectory(dep/eigen)

# main project
add_subdirectory(src)