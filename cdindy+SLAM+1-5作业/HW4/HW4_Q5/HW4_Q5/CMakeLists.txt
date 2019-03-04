project(HW4_Q5)
cmake_minimum_required(VERSION 2.8)


set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-O3" )

#SET c++
set( CMAKE_CXX_FLAGS "-std=c++11" )

#find opencv3
find_package( OpenCV 3 REQUIRED )


#include
include_directories(
${OPENCV_INCLUDE_DIRS}   #opencv
${EIGEN3_INCLUDE_DIR}    #Eigen
)


add_executable( gaussnewton gaussnewton.cpp)

target_link_libraries( gaussnewton
${OpenCV_LIBS}
${EIGEN3_LIBS}
)

