# --------------------------------------------------------------------------- #
#    CMake file for UCBlock/tools                                             #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----- nc4generator -------------------------------------------------------- #
add_executable(nc4generator nc4generator.cpp)
target_compile_features(nc4generator PRIVATE cxx_std_17)
target_link_libraries(nc4generator PRIVATE SMS++::SMS++)

# ----- Install instructions ------------------------------------------------ #
include(GNUInstallDirs)
install(TARGETS nc4generator
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# --------------------------------------------------------------------------- #
