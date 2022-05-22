# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- UnitBlock can be scaled (replicated).
- BatteryUnitBlock, IntermittentUnitBlock, and ThermalUnitBlock implement scale().
- BatteryUnitBlock and IntermittentUnitBlock can have their minimum and
  maximum power and storage levels scaled.

### Changed

- Improved UCBlock abstract constraints code.

### Fixed

- Serialization of BatteryUnitBlock, HydroUnitBlock, IntermittentUnitBlock,
  NetworkBlock, ThermalUnitBlock, UCBlock.

- Deserialization of DCNetworkBlock and NetworkBlock.

## [0.6.0] - 2021-12-08

### Added

- ThermalUnitDPSolver.
- Option to add spinning reserve variables to the Objective of ThermalUnitBlock.
- is_feasible() to UnitBlocks.

### Fixed

- Bugs in deserialization.
- Bug in the Objective of ThermalUnitBlock.
- Bugs in Constraints.
- Bugs in methods that are used to retrieve data.
- Bugs in methods to change the physical representation.

## [0.5.0] - 2021-05-02

### Fixed

- Too many fixes to list.

## [0.4.1] - 2020-09-16

### Fixed

- Generation of abstract Constraint in UCBlock

## [0.4.0] - 2020-09-16

### Added

- Support for Hydro[System]UnitBlock.

## [0.3.1] - 2020-07-17

### Fixed

- Removed a bugged test file.

## [0.3.0] - 2020-07-15

### Added

- Methods for changing data after deserialization.
- ThermalUnitBlock unit tests.
- Added all the UCBlock codes.

### Changed

- Some getter methods for ThermalUnitBlock data.

## [0.2.0] - 2020-03-06

### Added

- HydroSystemUnitBlock.
- Conan recipe.

### Fixed

- Minor bugs.

## [0.1.0] - 2020-02-06

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/ucblock/-/compare/0.6.0...develop
[0.6.0]: https://gitlab.com/smspp/ucblock/-/compare/0.5.0...0.6.0
[0.5.0]: https://gitlab.com/smspp/ucblock/-/compare/0.4.1...0.5.0
[0.4.1]: https://gitlab.com/smspp/ucblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/ucblock/-/compare/0.3.1...0.4.0
[0.3.0]: https://gitlab.com/smspp/ucblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/ucblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/ucblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/ucblock/-/tags/0.1.0
