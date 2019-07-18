add_executable(nc4generator nc4generator.cpp helper.cpp)
target_link_libraries(nc4generator
        SMS++
        netcdf
        netcdf-cxx4)
