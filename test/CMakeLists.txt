find_package(netCDF REQUIRED)

add_executable(UCBlock_test ucblock_test.cpp)
target_include_directories(UCBlock_test PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(UCBlock_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(UCBlock_test
                      SMS++
                      UCBlock
                      netcdf
                      netcdf-cxx4
                      )

add_executable(ThermalUnitBlock_test thermalunitblock_test.cpp)
target_include_directories(ThermalUnitBlock_test PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(ThermalUnitBlock_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(ThermalUnitBlock_test
                      SMS++
                      UCBlock
                      MILPSolver
                      netcdf
                      netcdf-cxx4
                      )

# FIXME: Workaround for CPLEX
if (OVERRIDE_CPLEX_LIB)
    target_link_libraries(ThermalUnitBlock_test ${OVERRIDE_CPLEX_LIB})
else ()
    target_link_libraries(ThermalUnitBlock_test cplex-library)
endif ()
