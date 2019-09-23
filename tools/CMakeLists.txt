find_package(netCDF REQUIRED)

# ----- nc4generator -----------------------------------------------------------
add_executable(nc4generator nc4generator.cpp helper.cpp)
target_link_directories(nc4generator PRIVATE ${netCDF_LIB_DIR})

target_link_libraries(nc4generator
                      SMS++::SMSpp
                      SMS++::UCBlock
                      netcdf
                      netcdf-cxx4
                      )
