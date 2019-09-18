find_package(netCDF REQUIRED)
include_directories(${netCDF_INCLUDE_DIR})
link_directories(${netCDF_LIB_DIR})

add_executable(nc4generator nc4generator.cpp helper.cpp)
target_link_libraries(nc4generator
        SMS++
        netcdf
        netcdf-cxx4)
