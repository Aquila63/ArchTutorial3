cmake_minimum_required(VERSION 2.8)
project(ArchTutorial3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mrtm -pthread")

set(SOURCE_FILES sharing.cpp helper.cpp helper.h)
add_executable(ArchTutorial3 ${SOURCE_FILES})
