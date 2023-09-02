cmake_minimum_required(VERSION 3.1...3.23)
project(ga_spgmt_uu)

set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("main.cpp")