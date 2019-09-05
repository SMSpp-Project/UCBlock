add_library(
        UCBlock OBJECT
        src/UCBlock.cpp
        src/NetworkBlock.cpp
        src/DCNetworkBlock.cpp
        src/BusNetworkBlock.cpp
        src/UnitBlock.cpp
        #src/HeatBlock.cpp
        src/ThermalUnitBlock.cpp
        src/HydroUnitBlock.cpp
        src/BatteryUnitBlock.cpp
        src/IntermittentUnitBlock.cpp)

add_subdirectory(tools)
add_subdirectory(test)
