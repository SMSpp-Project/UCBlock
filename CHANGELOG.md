# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `ThermalUnitBlock::set_min_up_down_time()` sets the minimum up and down
  times before the Variable are generated, which is when they can still be
  set: unlike the other data they decide how many Variable there are, and the
  method throws once those exist.

- The PTDF matrix of a `DCNetworkBlock` is no longer perturbed by a Tikhonov
  term on the diagonal (`f_tikhonov_coeff` is 0 by default), the reduced
  Laplacian being nonsingular once a reference node per connected component is
  removed. The perturbation broke the conservation of the flows: the PTDF rows,
  the nodal balance of the nodes an HVDC line touches and the overall balance
  are then no longer consistent, and together they pin a relation among the
  injections that has nothing to do with the network, so that on a network of
  both kinds of lines the PTDF formulation gave an optimum far from that of the
  KIRCHHOFF one, and could even declare the problem infeasible. The two now
  agree, on a network of 8 countries as on a triangle.

- Two lines that join the same two nodes sum their susceptances in the PTDF
  matrix, instead of the second one overwriting the entry of the first.

- The costs of the operating rules of a `NuclearUnitBlock`, i.e., those of a
  downward modulation step and of a deep decrease, can change:
  `set_down_modulation_costs()` and `set_deep_decrease_costs()` do it, each in
  its Range and Subset form and both in the methods factory, and a
  change of the corresponding coefficients of the Objective is folded into
  them, so that the physical representation follows and the DP Solver hears of
  it. They used to be refused, which stopped any Solver that writes on the
  Objective of the unit, such as the primal proximal heuristic of
  `LagrangianDualSolver` with its penalty term.

- `ThermalUnitBlock` accepts a `MinUpTime` (`MinDownTime`) of one instant more
  than the horizon, which says that the unit, being on (off) before it, never
  switches within it: the commitment is then fixed at every instant and there
  is no start-up or shut-down variable. The bound was the horizon itself, so
  that the last instant was free however long the minimum time.

- `ThermalUnitBlock` reads the optional variable `ShutDownCost`, the cost the
  unit pays at the instant it goes off, mapped over the time horizon as
  `StartUpCost` is and changed with `set_shutdown_costs()`. A unit that pays
  nothing to shut down keeps the vector empty and its Objective has no term
  for the shut-down variables, as before; the two DP Solvers charge the cost
  on the arc that closes an ON run.

- the pollutant budget constraints of `UCBlock` can also bound the emission
  from below (`PollutantMinBudget`, and with an equal upper bound match a
  value, changed with `set_pollutant_min_budget()`) and take the levels of
  the storages into account (`PollutantStorageRho`, a factor on the level of
  each storage at each time, the storages of a unit being at the node of its
  first generator), which is how a limit that charges a storage for the
  change of its level is written. To this end `UnitBlock` has
  `get_number_storages()` and `get_storage_level()`, returning 0 and nullptr
  unless redefined, as `BatteryUnitBlock` (its charge), `HydroUnitBlock` (its
  reservoirs) and `HydroSystemUnitBlock` (the reservoirs of all its units) do

- the pollutant budget constraints of `UCBlock`: for each zone of each
  pollutant, the emission of the electrical generators at the nodes of the
  zone summed over the whole time horizon, i.e., the conversion factor
  `PollutantRho` (which may depend on time, and includes the duration of the
  time instant) times the active power, cannot exceed `PollutantBudget`. The
  constraints are one vector with the zones of all the pollutants one after
  the other, as the budgets are in the file, they follow the scaling of the
  units, their duals are in the `UCBlockSolution` (bit 128), and the budget
  can be changed with `set_pollutant_budget()`; if `NumberPollutantZones` is
  not given every pollutant has one zone, and if `PollutantZones` is not given
  then all the nodes are in it. The constraints are grouped by pollutant, each
  with as many as its zones (`get_pollutant_constraints()[ p ][ z ]`), while
  the budgets and the duals keep the zones of all the pollutants one after
  the other, as the file does

### Fixed

- `NuclearUnitBlock::set_solution()` derived the start of a modulation and the
  three indicators of a deep decrease, and left the end of a modulation and
  the band of the output as they were, so that a Solver that writes the
  canonical Variable alone left them with the values of whoever had written
  them before: the schedule that came out then broke `ModulationEnd` by one
  and, since the band rows are relaxed by `MaxPower - MinPower` while a
  modulation travels, `BandPower` by that much. Both are derived now, the
  band by a sweep over the three of them, since it only changes where a
  modulation ends and a band that holds the output at one instant may leave
  none at the next

- The band of a `NuclearUnitBlock` could move against the direction of the
  modulation that moves it, the rows asking only that it change by one: where
  the landing power is a breakpoint, and both bands hold it, a modulation
  upwards could therefore leave the unit a band lower, which costs nothing
  while a downward one is paid for. The move now follows the direction, as it
  does in the dynamic program, the two having disagreed by up to `1.2e-4` on
  the instances of the battery. `has_modulation_direction()` is true whenever
  the output is banded, the direction Variable being what says which way the
  band goes

- At the last instant of the horizon a `NuclearUnitBlock` with the bands could
  take a modulation step that was neither a full ramp nor landed in a band:
  the full ramp of a step is imposed through the modulation of the next
  instant, which is not there, while the end of the modulation, which the band
  rows are relaxed without, was free. The step was therefore a last one for the
  ramp and not for the band, and on one instance of the battery the unit kept a
  lower output for the whole day and then modulated by a fraction of the ramp
  at its end, `5e-6` below what the dynamic program, where a modulation either
  ends in a band or goes on by the full ramp, could reach. At the last instant
  the next modulation is now replaced by the one that goes on past the horizon,
  i.e., `m_{T-1} - e_{T-1}`, a single-step modulation always ends there, and
  `set_solution()` marks as ended a modulation whose last step is not a full
  ramp

- `ACNetworkBlock` registered `"DCNetworkBlock::set_active_demand"` again,
  with an adapter casting to `ACNetworkBlock`, so that which of the two
  registrations survived depended on the order of the static initialization,
  and a call by name on a `DCNetworkBlock` could go through a cast to a class
  it is not. `DCNetworkBlock` registers it, and it reaches an
  `ACNetworkBlock` as well

- `DCNetworkBlock` scaled its flow limits by `C_v_scal` only in the box of
  the lines without a design variable, and lost the factor there too as soon
  as `set_kappa()` rewrote them: the rows of the lines with a design variable
  now carry it as well, and keep it through `set_kappa()`. The default factor
  is 1, which is why nothing showed

- `DCNetworkBlock` gave a line with a design variable its lower row only when
  kappa times the minimum flow was not zero, so that a line starting with a
  zero kappa and a nonzero minimum flow was left with the lower half of the
  box, i.e., with a nonnegative flow, and `set_kappa()` had no row to write
  the new lower limit in: the row now depends on the minimum flow alone, as
  the documentation of `v_design_min_row` already said

- `BatteryUnitBlock::get_max_power_constraints()` and
  `get_max_outtake_binary_constraints()` returned the second element of the
  first row of their two-row group rather than the first element of the
  second, so that `get_kappa_linearization()` read, for every instant, the
  multiplier of the row of the instant after it, and the coefficient the
  InvestmentFunction gets for a battery under investment was that of the
  wrong rows. The coefficient is now also read off one side at a time: where
  kappa is 0 the two rows pin the variable and both multipliers can be
  nonzero on what is a single equality, and summing the two terms made the
  linearization steeper than the function it linearizes

- `HydroUnitBlock` gave an arc that is neither a turbine nor a pump at some
  instant (no flow allowed) a single flow-to-power row, but then went on as if
  it had one per piece: with more than one piece on that arc, the rows of the
  arcs after it at that instant took the linear and constant terms of other
  pieces. `FlowActivePower_Const` is now grouped by instant and arc, the rows
  of an arc being its own

- `DCNetworkBlock` had a flow definition row for every line, and left those
  of the HVDC lines without a Function, which `is_feasible()` cannot compute;
  `ACNetworkBlock` had the two product-of-voltages variables for every line,
  and used them only on the DC lines; `DesignNetworkBlock` had a bound for
  every design variable, and left empty those of the variables fixed to 1.
  They now have them only where they are used

- `UCBlock::is_feasible()` answered false on an optimal solution: it checked
  the sub-Block with no Configuration, i.e., with tolerance 0 whatever the
  tolerance of the UCBlock, and `DCNetworkBlock::is_feasible()` checked the
  overall balance constraint also in the formulations that do not generate
  it, where a constraint with no Function cannot be computed. The sub-Block
  of `UCBlock`, `DesignNetworkBlock` and `HydroSystemUnitBlock` are now
  checked with the tolerance and type of violation of their father, unless
  their BlockConfig has a Configuration of its own, and the overall balance
  only where it exists

- `SecondaryZones`, `InertiaZones` and `NumberPollutantZones` were only read
  if the demand of the previous family of constraints was in the file

- `UCBlock::add_Modification` looked for the scaling of a unit only at the
  first level of the Modification it is given: the scaling of several units
  can arrive inside one GroupModification, and there the rows that use the
  Variable of those units, i.e., the node injection and the demand ones,
  were left untouched

### Changed

- `NuclearUnitBlock` keeps each family of operating rules in a group of its
  own (the tight rows on the starts of a modulation, `ModulationStartsApart`,
  `ModulationEndStarts` and `ModulationStepStarted`, are no longer mixed with
  `ModulationStability` and `ModulationMaxLength`), and the groups whose
  number of rows depends on the instant (`StartUpStability`, `BandKeep`,
  `BandMove`, `ModulationEndLink`) have the rows of instant t in their entry
  t

- the setters that change one datum spanning the whole time horizon issue
  their "abstract" Modification inside a GroupModification, one per setter,
  rather than one loose Modification per instant: a Solver able to execute a
  whole set of changes in one operation can then do so, while one that is not
  takes the group apart and sees exactly what it saw before. So far
  `SlackUnitBlock::set_active_power_cost`, `HydroUnitBlock::set_inflow`,
  `HydroUnitBlock::update_initial_flow_rate_in_cnstrs`,
  `ThermalUnitBlock::set_maximum_power`,
  `IntermittentUnitBlock::update_max_power_in_cnstrs` and `set_kappa`,
  `BatteryUnitBlock::update_kappa_in_cnstrs`,
  `DCNetworkBlock::change_power_flow_limit_constraints`, the four setters of
  the prices and of the demand of `ECNetworkBlock`, and the reaction of
  `UCBlock` to the scaling of a unit, where the four `update_*_constraints`
  now travel in one channel, `DCNetworkBlock::set_network_cost`, the
  `set_active_power_cost` of `HydroUnitBlock` and of
  `IntermittentUnitBlock`, `IntermittentUnitBlock::update_objective` and
  `ThermalUnitBlock::update_objective_active_power`. The last two were
  `const`, which opening a channel is not: they are private helpers called
  only from methods that are not const, and the const is gone

## [0.9.0] - 2026-09-13

### Added

- the operating rules of a nuclear unit in the `NuclearUnitBlock`:
  modulations of up to `MaxModulationLength` steps at the full ramp with the
  last one free, the stability of `ModulationTime` instants that follows a
  modulation and the `StabilityAfterStartUp` instants that follow a start-up,
  the daily limits (`DayLength`, `ModulationsPerDay`, `StartUpsPerDay`,
  `DeepDecreasesPerDay`), the deep decreases (`DeepDecreaseThreshold`,
  `DeepDecreaseGradient`, `DeepDecreaseCost`), the cost of a downward
  modulation step (`DownModulationCost`) and the splitting of the output
  into `PowerBands`, each written only when the data asks for it, so that a
  file of the previous model keeps its meaning

- a tight version of the rules, selected by two further bits of the same
  integer `Configuration` that chooses the formulation (`TightRules`,
  `TightRamp` and `TightCuts`), the third of which gives the tight rows by
  separation, in `generate_dynamic_constraints()`, rather than writing them

- `NuclearUnitExtDPSolver` rewritten as a labelled run-length dynamic
  program over the `ThermalUnitExtDPSolver`, which grows the hooks the
  labels need (`on_moves()`, `shut_label()`, `idle_label()`,
  `start_labels()`, `label_dominates()`, `trim_domination()`)

## [0.8.0] - 2026-09-12

### Added

- `get_Solution()` for the three dynamic programming Solver of the
  ThermalUnitBlock, which fill the ThermalUnitBlockSolution straight out of
  the schedule the DP has found, without writing anything into the Block and
  therefore without requiring any Variable to exist; the recovery of the
  schedule, which get_var_solution() shares, is factored into
  `recover_schedule()`

- accessors and setters to the parts of the solution saved in
  `UnitBlockSolution` and to the dimensioning variable saved in
  `ThermalUnitBlockSolution`

- the factor the flow limits of a DCNetworkBlock are scaled by is readable,
  since whoever reads their duals needs it

### Changed

- an auxiliary Variable, and the rows that fence it, are only generated when
  the data asks for them: the intake and the outtake level of the
  `BatteryUnitBlock` when the storing and the extracting efficiency disagree
  at some time instant, the linearisation of the flow cost of the
  `DCNetworkBlock` and of the `OTSNetworkBlock` when some line is priced, the
  minimum power row of the `IntermittentUnitBlock` and the maximum and minimum
  power rows of the `HydroUnitBlock` when the unit produces reserve. Where the
  rows are not generated, what they reduce to is a bound on the active power,
  and it is stated as a bound

- the flow-to-power relation of a `HydroUnitBlock` arc is an equality, rather
  than the concave outer approximation a piecewise arc needs, when the arc has
  a single piece with no constant term and the reservoir it leaves has a
  spillage outlet: the flow can then be substituted away

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- `UnitBlockSolution::write()` threw on a generator having no Variable for a
  part the Solution carries, while `read()` skips it: since what a Solution
  is made of is now decided by the data, that part can be there while the
  Variable are not, and writing on what is there is the only reading of it
  consistent with the other direction

- `HydroSystemUnitBlock` said it has the primary and the secondary reserve
  whenever the enclosing UCBlock asks for them, whatever the HydroUnitBlock
  it is made of have: it now answers for the units, which need not agree

- the dynamic programming Solver of the ThermalUnitBlock crashed on a unit
  whose power domain collapses to a single point at some time instant, say
  min_power == max_power: the sliding minimum dropped the zero-width piece,
  the value function came out empty, and the caller read that as no piece at
  all. The point is now emitted as a zero-width piece, which the rest of the
  machinery already handles, and a transition that really is infeasible
  closes the run instead of being walked into

- the documentation of the formulations of the ThermalUnitBlock, which
  numbered the p_t one 3 and the dynamic programming one 2, the other way
  round with respect to ptForm and DPForm and to what set_variables() does

- the intake and outtake bounds of a BatteryUnitBlock were held in a
  LB0Constraint, whose LHS is fixed at zero by the type: right for the pair
  of one-sided fences on the intake and the outtake, wrong once the pair is
  folded onto the signed active power, where the bound is two-sided and its
  charging side is negative. Generating the constraints of a battery in the
  folded form threw "cannot change LHS in a LB0Constraint"; they are now a
  BoxConstraint, which the split form uses with its default LHS of zero

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

- `tools/csv2netCDF` from Energy Community Julia codebase

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

[Unreleased]: https://gitlab.com/smspp/ucblock/-/compare/0.8.0...develop
[0.8.0]: https://gitlab.com/smspp/ucblock/-/compare/0.7.0...0.8.0
[0.7.0]: https://gitlab.com/smspp/ucblock/-/compare/0.6.3...0.7.0
[0.6.3]: https://gitlab.com/smspp/ucblock/-/compare/0.6.2...0.6.3
[0.6.2]: https://gitlab.com/smspp/ucblock/-/compare/0.6.1...0.6.2
[0.6.1]: https://gitlab.com/smspp/ucblock/-/compare/0.6.0...0.6.1
[0.6.0]: https://gitlab.com/smspp/ucblock/-/compare/0.5.0...0.6.0
[0.5.0]: https://gitlab.com/smspp/ucblock/-/compare/0.4.1...0.5.0
[0.4.1]: https://gitlab.com/smspp/ucblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/ucblock/-/compare/0.3.1...0.4.0
[0.3.1]: https://gitlab.com/smspp/ucblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/ucblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/ucblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/ucblock/-/tags/0.1.0
