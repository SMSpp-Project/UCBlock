# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ UCBlock                                        #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
project(UCBlock VERSION 0.4.1 LANGUAGES CXX)

# ----- Requirements -------------------------------------------------------- #
find_package(SMSpp REQUIRED)

# ----- Settings ------------------------------------------------------------ #
# Debug checks and printouts
if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_compile_definitions(NDEBUG)
endif ()

# ----- Library ------------------------------------------------------------- #
add_library(UCBlock )
target_compile_features(UCBlock PUBLIC cxx_std_17)

target_sources(UCBlock PRIVATE
               src/UCBlock.cpp
               src/NetworkBlock.cpp
               src/DCNetworkBlock.cpp
               src/BusNetworkBlock.cpp
               src/UnitBlock.cpp
               #src/HeatBlock.cpp
               src/BatteryUnitBlock.cpp
               src/ThermalUnitBlock.cpp
               src/HydroUnitBlock.cpp
               src/IntermittentUnitBlock.cpp
               src/SlackUnitBlock.cpp
               src/HydroSystemUnitBlock.cpp
               )

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        UCBlock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(UCBlock PUBLIC
                      SMS++::SMSpp)

# Alias for testing
add_library(SMS++::UCBlock ALIAS UCBlock)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS UCBlock
        EXPORT UCBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT UCBlock-targets
        FILE UCBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/UCBlock)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/UCBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/UCBlock)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/UCBlock)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    export(EXPORT UCBlock-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/UCBlockTargets.cmake
           NAMESPACE SMS++::)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE UCBlock)
endif ()

# ----- Tools --------------------------------------------------------------- #
add_subdirectory(tools)

# ----- Unit testing -------------------------------------------------------- #
#if (BUILD_TESTING)
#    add_subdirectory(test)
#endif()
