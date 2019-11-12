# --------------------------------------------------------------------------- #
#    Main CMake file for UCBlock tools                                        #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
# TODO: set proper version
project(UCBlock_tools VERSION 1.0.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
# Since we are using the block factory, objects from statically linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.

if (NOT BUILD_SHARED_LIBS)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        add_link_options("-Wl,-all_load")
    else()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

# ----- Requirements -------------------------------------------------------- #
find_package(MILPSolver REQUIRED)
find_package(UCBlock REQUIRED)

# ----- nc4generator -------------------------------------------------------- #
add_executable(nc4generator nc4generator.cpp)
target_link_libraries(nc4generator PRIVATE SMS++::UCBlock)

# ----- thermalunit_solver -------------------------------------------------- #
add_executable(thermalunit_solver thermalunit_solver.cpp)
target_link_libraries(thermalunit_solver PRIVATE
                      SMS++::UCBlock
                      SMS++::MILPSolver)

# ----- ucblock_solver ------------------------------------------------------ #
add_executable(ucblock_solver ucblock_solver.cpp)
target_link_libraries(ucblock_solver PRIVATE
                      SMS++::UCBlock
                      SMS++::MILPSolver)
