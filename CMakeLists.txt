add_library(
        UCBlock OBJECT
        src/UCBlock.cpp
        src/NetworkBlock.cpp
        src/DCNetworkBlock.cpp
        src/BusNetworkBlock.cpp
        src/UnitBlock.cpp
        src/HeatBlock.cpp
        src/ThermalUnitBlock.cpp
        src/HydroUnitBlock.cpp
        src/BatteryUnitBlock.cpp
        src/CentralizedDemandResponseUnitBlock.cpp
        src/EMobilityUnitBlock.cpp
        src/IntermittentGenerationUnitBlock.cpp
        src/PowerToGasUnitBlock.cpp)

add_subdirectory(tools)
