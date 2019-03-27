project(HW7_Q3)
cmake_minimum_required(VERSION 2.8)

set( CMAKE_BUILD_TYPE "debug" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

# 寻找G2O Cholmod eigen3
find_package( G2O REQUIRED )
find_package( Cholmod )
include_directories(
    ${G2O_INCLUDE_DIRS} ${CHOLMOD_INCLUDE_DIR}
    "/usr/include/eigen3"
)

# sophus
find_package( Sophus REQUIRED )
include_directories( ${Sophus_INCLUDE_DIRS} )

#OpenCV
find_package( OpenCV 3 REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

# find Pangolin
find_package(Pangolin REQUIRED)
include_directories( ${Pangolin_INCLUDE_DIRS} )



SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear
    g2o_types_icp g2o_types_slam2d g2o_types_sba g2o_types_slam3d g2o_core g2o_interface
    g2o_solver_csparse g2o_solver_structure_only g2o_csparse_extension g2o_opengl_helper g2o_solver_dense
    g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )


add_executable(${PROJECT_NAME} directBA.cpp)

target_link_libraries( ${PROJECT_NAME}

        ${EIGEN3_LIBS}
        ${G2O_LIBS} ${CHOLMOD_LIBRARIES}
        ${CHOLMOD_LIBRARIES}
        ${Sophus_LIBRARIES}
        ${OpenCV_LIBRARIES}
        ${Pangolin_LIBRARIES}
        )
