find_package(netCDF REQUIRED)
include_directories(${netCDF_INCLUDE_DIR})
link_directories(${netCDF_LIB_DIR})

add_executable(UCBlock_test ucblock_test.cpp)
target_link_libraries(UCBlock_test
        SMS++
        UCBlock
        netcdf
        netcdf-cxx4)

add_executable(ThermalUnitBlock_test thermalunitblock_test.cpp)
target_link_libraries(ThermalUnitBlock_test
        SMS++
        UCBlock
        MILPSolver
        cplex-library
        netcdf
        netcdf-cxx4)
