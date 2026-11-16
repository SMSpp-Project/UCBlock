# --------------------------------------------------------------------------- #
#    CMake file for UCBlock/test                                              #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# ----- UCBlock_unit_test --------------------------------------------------- #
add_executable(${modName}_unit_test test.cpp test_power_limits.cpp)

# test.cpp and test_power_limits.cpp read their configuration files (the
# BlockConfig of the formulations, TUBCfg-*.txt, NUBCfg-*.txt, ACBCfg-*.txt
# and PFBCfg.txt, the BlockSolverConfig LPCmp*.txt and the Configuration of
# the Solution *SolCfg-*.txt) from the source directory
target_compile_definitions(${modName}_unit_test PRIVATE
                           UCBLOCK_TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${modName}_unit_test PRIVATE SMS++::UCBlock)

# The signs of the dual values are checked with a :MILPSolver when the build
# has one; MILPSolver comes after this module in the umbrella, hence the
# option is looked at as well as the target. The module alone has none, and
# the tester then skips that check.
if (TARGET SMS++::MILPSolver OR BUILD_MILPSolver)
    target_link_libraries(${modName}_unit_test PRIVATE SMS++::MILPSolver)
endif ()

smspp_copy_runtime(${modName}_unit_test)

# The tester builds its instances in memory and needs no data; it writes its
# netCDF groups in a file of the build directory and removes it.
add_test(NAME ${modName}_unit_test
         COMMAND ${modName}_unit_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ----- the label of this module's tests ------------------------------------ #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}")
endif ()

# --------------------------------------------------------------------------- #
