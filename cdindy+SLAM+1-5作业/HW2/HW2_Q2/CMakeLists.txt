project(HW2)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-O3" )

include_directories( "/usr/include/eigen3" )

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

