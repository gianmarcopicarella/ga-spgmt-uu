set(CGAL_DATA_DIR ${PROJECT_SOURCE_DIR}/data)

find_package(CGAL REQUIRED)

file(   GLOB_RECURSE 
        COMMON
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp)

file(   GLOB_RECURSE 
        MAIN_ALL
        ${CMAKE_CURRENT_SOURCE_DIR}/main/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/main/*.cpp)

file(   GLOB_RECURSE 
        TESTS_ALL
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)

add_executable(main ${MAIN_ALL} ${COMMON})
add_executable(tests ${TESTS_ALL} ${COMMON})

if ( CGAL_AUTO_LINK_ENABLED )
target_link_libraries(main ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(tests ${CGAL_3RD_PARTY_LIBRARIES} )
else()
target_link_libraries(main ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(tests ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md
target_link_libraries(main PRIVATE CGAL::CGAL CGAL::Data)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain CGAL::CGAL CGAL::Data)