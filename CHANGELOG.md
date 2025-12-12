# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

### Changed 

### Fixed 

## [0.7.0] - 2025-12-12

### Added 

- DesignNetworkBlockSolution

- [big] DesignNetworkBlock allowing to do desing of lines across the
  time horizon

- simple stochastic case in csv2netcdf

- constant term to UCBlock

- design variables in Solution for BatteryUnitBlock,
  IntermittentUnitBlock, ThermalUnitBlock

- MinCapacityDesign for IntermittentUnitBlock and BatteryUnitBlock

- capacity expansion for Intermittent and Battery units

- [big] "compact" NetworkBlockSolution to avoid issue with netCDF
  files being very slow when writing many small sub-group

- [huge] hyperarcs in DCNetworkBlock, our quick-and-dirty solution
  to multi-energy modelling

- [big] efficiency of lines in DCNetworkBlock

- HydroUnitBlockSolution, BatteryUnitBlockSolution,
  HydroSystemUnitBlockSolution

- cycling for hydro (volume at the end of the interval equal to that
  at the beginning)

- add ActivePowerCost to HydroUnitBlock

- Include LinearTerm to IntermittentUnitBlock

- DCNetworkBlockSolution

- UCBlockSolution, UnitBlockSolution, NetworkBlockSolution

### Changed 

- type of ActivePower\_Bound\_Const for slack unit

- allow negative IntermittentUB and negative loads

- major data handling upgrade whereby instances are no longer included
  in the repo but can be downloaded

- refactored messy HydroUnitBlock constraints building

- checked the value of delta ramp-up/down

- moved node injection bound constraints to NetworkBlock

### Fixed 

- UCBlock::set\_active_power\_demand methods

- separated battery desing from converter design in BatteryUnitBlock

- IntakeOuttake\_Design\_Battery constraint in BatteryUnitBlock

- Demand Battery constraints in BatteryUnitBlock

- issue in secondary and inertia demand constraints

- get\_min\_power() and get\_max\_power-89 in HydroUnitBlock

- bug in DCNetworkBlock::change\_power\_flow\_limit\_constraints()

- catastrophic blunder in NuclearUnitBlock

- bugs when getting AC/DC lines

- issue in HydroUnitBlock when having 0 pieces

- is\_cooperative in ECNB

- cleaned up factory issue in NetworkBlock

- objective in ECNB

- many minor fixes

## [0.6.3] - 2024-02-29

### Added 

- design variables in ThermalUnitBlock, BatteryUnitBlock,
  IntermittentUnitBlock

- tools/DataConverter from Energy Community Julia codebase

- `data/nc4/EC_Data` test data sets

- ECNetworkBlock

### Changed 

- adapted to new CMake / makefile organisation

- NetworkBlock can now span multiple time instants

- updated Julia and nc4 files with the stochastic logic

### Fixed

- bug in `ThermalUnitBlock::update_objective_start_up()` in which
  the `v_start_up` vector was being accessed at a wrong index

- bugs when retrieving and checking constraints in BatteryUnitBlock

- bug in ThermalUnitBlock::update\_objective\_start\_up()

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

### Fixed

- minor bugs

## [0.1.0] - 2020-02-06

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/ucblock/-/compare/0.7.0...develop
[0.7.0]: https://gitlab.com/smspp/ucblock/-/compare/0.6.3...0.7.0
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
