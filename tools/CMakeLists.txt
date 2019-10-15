find_package(netCDF REQUIRED)
find_package(CPLEX)

# ----- nc4generator --------------------------------------------------------- #
add_executable(nc4generator nc4generator.cpp)

target_link_libraries(nc4generator PRIVATE
                      SMS++::UCBlock
                      NetCDF::NetCDFCxx)

if (CPLEX_FOUND)
    add_executable(thermalunit_solver thermalunit_solver.cpp)
    target_link_libraries(thermalunit_solver PRIVATE
                          SMS++::UCBlock
                          SMS++::MILPSolver
                          CPLEX::Cplex)
endif ()
