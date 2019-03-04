project(HW3)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-O3" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3 -march=native")

# Check C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

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
