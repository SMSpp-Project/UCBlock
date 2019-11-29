# ----- nc4generator -------------------------------------------------------- #
add_executable(nc4generator nc4generator.cpp)
target_link_libraries(nc4generator PRIVATE SMS++::UCBlock)
