# ----- Project description and directories ------------------------------------
cmake_minimum_required(VERSION 3.8)
project(UCBlock
        VERSION 1.0.0 # TODO: set proper version
        HOMEPAGE_URL "https://gitlab.com/AliGhezel/ucblock"
        DESCRIPTION "UC Blocks for SMS++")

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SMSpp REQUIRED)

# ----- Library ----------------------------------------------------------------
add_library(UCBlock STATIC)
target_compile_features(UCBlock PUBLIC cxx_std_17)

target_sources(UCBlock PRIVATE
    ${PROJECT_SOURCE_DIR}/UCBlock.cpp
    ${PROJECT_SOURCE_DIR}/NetworkBlock.cpp
    ${PROJECT_SOURCE_DIR}/DCNetworkBlock.cpp
    ${PROJECT_SOURCE_DIR}/BusNetworkBlock.cpp
    ${PROJECT_SOURCE_DIR}/UnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/HeatBlock.cpp
    ${PROJECT_SOURCE_DIR}/ThermalUnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/HydroUnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/BatteryUnitBlock.cpp
    ${PROJECT_SOURCE_DIR}/IntermittentUnitBlock.cpp
    )

target_include_directories(UCBlock PUBLIC
                           $<INSTALL_INTERFACE:include/SMS++>
                           $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>)

target_link_libraries(UCBlock PUBLIC
                      SMS++::SMSpp
                      Boost::boost
                      Eigen3::Eigen)

# Alias for testing
add_library(SMS++::UCBlock ALIAS UCBlock)

# ----- Subdirectories ---------------------------------------------------------
#add_subdirectory(tools)
#add_subdirectory(test)

# ----- Install instructions ---------------------------------------------------
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Install library
install(TARGETS UCBlock
        EXPORT  UCBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install export file *Targets.cmake
install(EXPORT UCBlock-targets
        FILE   UCBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${INSTALL_CONFIGDIR})

# Create *ConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Create *Config.cmake
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/UCBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

# Install the Config and ConfigVersion
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/UCBlockConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR})

# ----- Export from the build tree ---------------------------------------------
export(EXPORT UCBlock-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/UCBlockTargets.cmake
       NAMESPACE SMS++::)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE UCBlock)

