cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

project(UCBlock DESCRIPTION "UC Blocks for SMS++")

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(UCBlock_SOURCES
    ${PROJECT_SOURCE_DIR}/UCBlock.cpp
    ${PROJECT_SOURCE_DIR}/NetworkBlock.cpp
    ${PROJECT_SOURCE_DIR}/DCNetworkBlock.cpp
    ${PROJECT_SOURCE_DIR}/BusNetworkBlock.cpp
    ${PROJECT_SOURCE_DIR}/UnitBlock.cpp
    #${PROJECT_SOURCE_DIR}/HeatBlock.cpp
    ${PROJECT_SOURCE_DIR}/ThermalUnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/HydroUnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/BatteryUnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/IntermittentUnitBlock.cpp
    )

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# TODO: Write and use Find*.cmake files
include_directories(../SMS++/include)

include_directories("${PROJECT_INCLUDE_DIR}")
add_library(${PROJECT_NAME} STATIC ${UCBlock_SOURCES})

add_subdirectory(tools)
add_subdirectory(test)
