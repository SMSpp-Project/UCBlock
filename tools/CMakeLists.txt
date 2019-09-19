find_package(netCDF REQUIRED)

add_executable(nc4generator nc4generator.cpp helper.cpp)
target_include_directories(nc4generator PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(nc4generator PRIVATE ${netCDF_LIB_DIR})

target_link_libraries(nc4generator
                      SMS++
                      UCBlock
                      netcdf
                      netcdf-cxx4
                      )
