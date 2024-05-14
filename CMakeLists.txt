cmake_minimum_required(VERSION 3.25)
project(DiagramDynamo)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
add_executable(DiagramDynamo main.cpp TransferFunctionBlock.cpp TransferFunctionBlock.h)
