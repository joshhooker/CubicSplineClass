cmake_minimum_required(VERSION 2.8)
project(spline)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
#set(CMAKE_VERBOSE_MAKEFILE ON)

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

add_executable(spline_test ${sources} ${headers})
