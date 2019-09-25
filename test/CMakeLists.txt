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
