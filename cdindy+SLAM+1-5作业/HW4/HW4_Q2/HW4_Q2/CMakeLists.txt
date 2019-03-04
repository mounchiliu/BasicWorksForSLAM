project(HW4)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-O3" )

#SET c++
set( CMAKE_CXX_FLAGS "-std=c++11" )

#find opencv 3
find_package( OpenCV 3 REQUIRED )
#include
include_directories( ${OPENCV_INCLUDE_DIRS} )



add_executable( undistort undistort_image.cpp )

target_link_libraries( undistort ${OpenCV_LIBS} )
