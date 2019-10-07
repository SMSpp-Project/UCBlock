## ----- Fetch GoogleTest ---------------------------------------------------- #
#include(FetchContent)
#fetchcontent_declare(
#        googletest
#        GIT_REPOSITORY https://github.com/google/googletest.git
#        GIT_TAG release-1.8.1
#)
#
#fetchcontent_getproperties(googletest)
#if (NOT googletest_POPULATED)
#    fetchcontent_populate(googletest)
#    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
#endif ()
#
## ----- ThermalUnitBlock_test ----------------------------------------------- #
#find_package(CPLEX)
#if (CPLEX_FOUND)
#
#    add_executable(ThermalUnitBlock_test tests_ThermalUnitBlock.cpp)
#
#    target_link_libraries(ThermalUnitBlock_test PRIVATE
#                          SMS++::UCBlock
#                          SMS++::MILPSolver
#                          NetCDF::NetCDFCxx
#                          CPLEX::Cplex
#                          gtest)
#
#    add_test(
#            NAME ThermalUnitBlock_test
#            COMMAND ThermalUnitBlock_test
#            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    )
#endif ()



# ----- UCBlock_test -------------------------------------------------------- #
add_executable(UCBlock_test ucblock_test.cpp)
target_link_libraries(UCBlock_test PRIVATE SMS++::UCBlock)

# ----- ThermalUnitBlock_test ----------------------------------------------- #
find_package(CPLEX)
if (CPLEX_FOUND)
    add_executable(ThermalUnitBlock_test thermalunitblock_test.cpp)
    target_link_libraries(ThermalUnitBlock_test PRIVATE
                          SMS++::UCBlock
                          SMS++::MILPSolver
                          CPLEX::Cplex)
endif ()
