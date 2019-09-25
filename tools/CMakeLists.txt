find_package(netCDF REQUIRED)

# ----- nc4generator --------------------------------------------------------- #
add_executable(nc4generator nc4generator.cpp helper.cpp)

target_link_libraries(nc4generator PRIVATE
                      SMS++::UCBlock
                      NetCDF::NetCDFCxx)
