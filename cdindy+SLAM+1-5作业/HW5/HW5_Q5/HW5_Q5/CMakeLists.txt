project(HW5_Q5)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_BUILD_TYPE "debug" )
set( CMAKE_CXX_FLAGS "-O3" )

set( CMAKE_CXX_FLAGS "-std=c++11" )

# find sophus
find_package( Sophus REQUIRED )
# find Pangolin
find_package(Pangolin REQUIRED)

include_directories(
${EIGEN3_INCLUDE_DIR}    #Eigen
${Sophus_INCLUDE_DIRS}   #Sophus
${Pangolin_INCLUDE_DIRS} #Pangolin
)

add_executable(${PROJECT_NAME} draw_trajectory.cpp)
target_link_libraries(${PROJECT_NAME}
${EIGEN3_LIBS}
${Pangolin_LIBRARIES}
${Sophus_LIBRARIES}
)
