project(HW6_Q2)
cmake_minimum_required(VERSION 2.8)
set( CMAKE_BUILD_TYPE "debug" )
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


add_executable( ${PROJECT_NAME} optical_flow.cpp)

target_link_libraries( ${PROJECT_NAME}
${OpenCV_LIBS}
${EIGEN3_LIBS}
)

