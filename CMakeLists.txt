add_library(
        UCBlock OBJECT
        src/UCBlock.cpp
        src/NetworkBlock.cpp
        src/DCNetworkBlock.cpp
        src/BusNetworkBlock.cpp
        src/UnitBlock.cpp
        src/ThermalUnitBlock.cpp
        src/HydroUnitBlock.cpp
        src/HeatBlock.cpp)

add_subdirectory(tools)
