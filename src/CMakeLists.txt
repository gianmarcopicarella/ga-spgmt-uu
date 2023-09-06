cmake_minimum_required(VERSION 3.1...3.23)
project(ga_spgmt_uu)

set(CMAKE_CXX_STANDARD 17)

set(CGAL_DATA_DIR "../data")

find_package(CGAL REQUIRED)

file(   GLOB_RECURSE 
        SRC_INCLUDE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

file(   GLOB_RECURSE 
        SRC_IMPL
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(main ${SRC_INCLUDE} ${SRC_IMPL})

if ( CGAL_AUTO_LINK_ENABLED )
target_link_libraries(main ${CGAL_3RD_PARTY_LIBRARIES} )
else()
target_link_libraries(main ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

target_link_libraries(main PUBLIC CGAL::CGAL CGAL::Data)

# Not working with multiple sources
#create_single_source_cgal_program("main.cpp")
#target_link_libraries(main PUBLIC CGAL::Data)