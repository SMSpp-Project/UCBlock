# ----- Fetch GoogleTest ---------------------------------------------------- #
include(FetchContent)
fetchcontent_declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.8.1
)

fetchcontent_getproperties(googletest)
if (NOT googletest_POPULATED)
    fetchcontent_populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif ()

# ----- Requirements -------------------------------------------------------- #
find_package(MILPSolver REQUIRED)
find_package(MCFBlock REQUIRED)
find_package(CPLEX)

# ----- ThermalUnitBlock_test ----------------------------------------------- #
if (CPLEX_FOUND)
    add_executable(ThermalUnitBlock_test tests_ThermalUnitBlock.cpp)
    target_link_libraries(ThermalUnitBlock_test PRIVATE
                          SMS++::UCBlock
                          SMS++::MILPSolver
                          CPLEX::Cplex
                          gtest)

    add_test(
            NAME ThermalUnitBlock_test
            COMMAND ThermalUnitBlock_test
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif ()
