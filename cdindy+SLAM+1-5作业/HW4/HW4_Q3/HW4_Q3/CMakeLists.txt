project(HW4_Q3)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-O3" )

#SET c++
set( CMAKE_CXX_FLAGS "-std=c++11" )

#find opencv3
find_package( OpenCV 3 REQUIRED )
# find Pangolin
find_package(Pangolin REQUIRED)
# find PCL
find_package(PCL REQUIRED)


#include
include_directories(
${OPENCV_INCLUDE_DIRS}   #opencv
${EIGEN3_INCLUDE_DIR}    #Eigen
${Pangolin_INCLUDE_DIRS} #Pangolin
${PCL_INCLUDE_DIRS}      #PCL
)


add_executable( disparity disparity.cpp )

target_link_libraries( disparity
${OpenCV_LIBS}
${EIGEN3_LIBS}
${Pangolin_LIBRARIES}
${PCL_LIBRARIES}
)
