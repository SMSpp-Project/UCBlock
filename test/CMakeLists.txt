find_package(netCDF REQUIRED)

# ----- UCBlock_test -----------------------------------------------------------
add_executable(UCBlock_test ucblock_test.cpp)

target_link_directories(UCBlock_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(UCBlock_test
                      SMS++::SMSpp
                      SMS++::UCBlock
                      netcdf
                      netcdf-cxx4
                      )

# ----- ThermalUnitBlock_test --------------------------------------------------
add_executable(ThermalUnitBlock_test thermalunitblock_test.cpp)

target_link_directories(ThermalUnitBlock_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(ThermalUnitBlock_test
                      SMS++::SMSpp
                      SMS++::UCBlock
                      SMS++::MILPSolver
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )
