add_executable(UCBlock_test ucblock_test.cpp)

target_link_libraries(UCBlock_test
        SMS++
        UCBlock
        netcdf
        netcdf-cxx4)
