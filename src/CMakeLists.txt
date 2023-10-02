set(CGAL_DATA_DIR ${PROJECT_SOURCE_DIR}/data)

find_package(CGAL REQUIRED)

find_package(Corrade REQUIRED Main)
find_package(Magnum REQUIRED GL
MeshTools
Primitives
Shaders
Trade
Sdl2Application)

# find_package (Eigen3 3.1.0 QUIET NO_MODULE)
# include(CGAL_Eigen3_support)

# IF IPE IS USED
# set(IPE_INSTALL_DIR "C:/Users/gi4nm/OneDrive/Desktop/ipe-7.2.28")
# set(IPE_INCLUDE_DIR ${IPE_INSTALL_DIR}/include)
# set(IPE_LIBRARIES ${IPE_INSTALL_DIR}/bin/ipe.dll)
# set(IPELET_INSTALL_DIR ${IPE_INSTALL_DIR}/ipelets)
# find_package(IPE 7)

# if(IPE_FOUND)
#     if(NOT ${IPE_VERSION} EQUAL "7")
#         message("NOTICE: ${IPE_VERSION} is not a supported version of IPE (only 7 is).")
#         set(IPE_FOUND FALSE)
#     endif()
# endif()

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

# file(   GLOB_RECURSE 
#         IPELETS_ALL
#         ${CMAKE_CURRENT_SOURCE_DIR}/ipelets/*.cpp)

add_executable(main ${MAIN_ALL} ${COMMON})
add_executable(tests ${TESTS_ALL} ${COMMON})

# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md
if ( CGAL_AUTO_LINK_ENABLED )
target_link_libraries(main ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(tests ${CGAL_3RD_PARTY_LIBRARIES} )
else()
target_link_libraries(main ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_link_libraries(tests ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

set(MAGNUM_MODULES Corrade::Main
Magnum::Application
Magnum::GL
Magnum::Magnum
Magnum::MeshTools
Magnum::Primitives
Magnum::Shaders
Magnum::Trade)

target_link_libraries(main PRIVATE CGAL::CGAL CGAL::Data ${MAGNUM_MODULES})
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain CGAL::CGAL CGAL::Data ${MAGNUM_MODULES})

# if(IPE_FOUND)
# foreach(IPELET_PATH ${IPELETS_ALL})
# get_filename_component(IPELET ${IPELET_PATH} NAME_WLE)
# add_library(CGAL_${IPELET} MODULE ${IPELET_PATH})
# target_include_directories(CGAL_${IPELET} BEFORE PRIVATE ${IPE_INCLUDE_DIR})
# target_link_libraries(CGAL_${IPELET} PRIVATE CGAL::CGAL CGAL::Eigen3_support ${IPE_LIBRARIES})
# if(IPELET_INSTALL_DIR)
# add_custom_command(TARGET CGAL_${IPELET}
#     POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:CGAL_${IPELET}> "${IPELET_INSTALL_DIR}"
#     COMMENT "Copy dll file to ${IPELET_INSTALL_DIR} directory" VERBATIM)
# endif ()
# endforeach(IPELET_PATH)
# endif()