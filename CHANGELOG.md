# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

### Changed 

### Fixed 

## [0.6.3] - 2024-02-29

### Added 

- design variables in ThermalUnitBlock, BatteryUnitBlock,
  IntermittentUnitBlock

- tools/DataConverter from Energy Community Julia codebase

- `netCDF_files/EC_Data` test data sets

- ECNetworkBlock

### Changed 

- adapted to new CMake / makefile organisation

- NetworkBlock can now span multiple time instants

- updated Julia and nc4 files with the stochastic logic

### Fixed

- bug in `ThermalUnitBlock::update_objective_start_up()` in which
  the `v_start_up` vector was being accessed at a wrong index

- bugs when retrieving and checking constraints in BatteryUnitBlock

- bug in ThermalUnitBlock::update_objective_start_up()

- separation of Perspective Cuts in ThermalUnitBlock

- the 3bin formulation including the start-up and shut-down limits
  cnstrs even if the ramp ones are not present

- default value for f_MinDownTime

- too many minor others to list

### Removed

- test/ moved to ThermalUnitBlock_Solver in test repository

- useless test_package


## [0.6.2] - 2023-05-17

### Added

- NuclearUnitBlock (didactic)

- is_feasible() to BatteryUnitBlock, SlackUnitBlock

- IntermittentUnitBlock::set_BlockConfig()

### Changed

- updated is_feasible() in HydroUnitBlock, DCNetworkBlock,
  IntermittentUnitBlock, ThermalUnitBlock, and BatteryUnitBlock

- removed "battery_type" from BatteryUnitBlock and check if negative prices may
  occur

## [0.6.1] - 2022-07-01

### Added

- UnitBlock can be scaled (replicated)

- BatteryUnitBlock, IntermittentUnitBlock, and ThermalUnitBlock implement
  scale()

- BatteryUnitBlock and IntermittentUnitBlock can have their minimum and
  maximum power and storage levels scaled (set_kappa())

### Changed

- improved UCBlock abstract constraints code

### Fixed

- serialization of BatteryUnitBlock, HydroUnitBlock, IntermittentUnitBlock,
  NetworkBlock, ThermalUnitBlock, UCBlock

- deserialization of DCNetworkBlock and NetworkBlock

## [0.6.0] - 2021-12-08

### Added

- ThermalUnitDPSolver

- Option to add spinning reserve variables to the Objective of ThermalUnitBlock

- is_feasible() to UnitBlocks

### Fixed

- bugs in deserialization

- bug in the Objective of ThermalUnitBlock

- bugs in Constraints

- bugs in methods that are used to retrieve data

- bugs in methods to change the physical representation

## [0.5.0] - 2021-05-02

### Fixed

- too many fixes to list

## [0.4.1] - 2020-09-16

### Fixed

- generation of abstract Constraint in UCBlock

## [0.4.0] - 2020-09-16

### Added

- support for Hydro[System]UnitBlock.

## [0.3.1] - 2020-07-17

### Removed

- a bugged test file

## [0.3.0] - 2020-07-15

### Added

- methods for changing data after deserialization

- ThermalUnitBlock unit tests

- all the UCBlock codes

### Changed

- some getter methods for ThermalUnitBlock data

## [0.2.0] - 2020-03-06

### Added

- HydroSystemUnitBlock

- Conan recipe

### Fixed

- minor bugs

## [0.1.0] - 2020-02-06

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/ucblock/-/compare/0.6.3...develop
[0.6.3]: https://gitlab.com/smspp/ucblock/-/compare/0.6.2...0.6.3
[0.6.2]: https://gitlab.com/smspp/ucblock/-/compare/0.6.1...0.6.2
[0.6.1]: https://gitlab.com/smspp/ucblock/-/compare/0.6.0...0.6.1
[0.6.0]: https://gitlab.com/smspp/ucblock/-/compare/0.5.0...0.6.0
[0.5.0]: https://gitlab.com/smspp/ucblock/-/compare/0.4.1...0.5.0
[0.4.1]: https://gitlab.com/smspp/ucblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/ucblock/-/compare/0.3.1...0.4.0
[0.3.0]: https://gitlab.com/smspp/ucblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/ucblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/ucblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/ucblock/-/tags/0.1.0
