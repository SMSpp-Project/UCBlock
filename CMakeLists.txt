cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(UCBlock DESCRIPTION "UC Blocks for SMS++")

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# ----- Custom includes --------------------------------------------------------
# TODO: When libraries will be installable, fix this
set(SMS++_INCLUDE_DIR ../SMS++/include)

# ----- Library ----------------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
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

add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${SMS++_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE_DIR})

# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(tools)
add_subdirectory(test)
