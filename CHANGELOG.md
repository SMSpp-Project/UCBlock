# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the unit test checks the continuous relaxations of the pt, SU, SD and SUSD
  formulations against the integer optimum on units on and off before the
  horizon, with both ramps, one or none, and the MILP of every formulation
  of a `NuclearUnitBlock`, with and without bands, against
  `NuclearUnitExtDPSolver`

- `NuclearUnitBlock` has seven more formulations of its operating rules,
  selected by bits 8-10 of the int `Configuration` that selects the
  formulation of the thermal part (`NuclearRules::FormMsk`): the one by the
  indicators of the states of the modulations (F0), the
  modulation-commitment one with the exact move of the output (F2) and its
  variant with the ends projected out (F3a), the one by the runs of the
  modulations (F4), and the flows on the label graph of the dynamic program
  without (F5) and with (F6) the counters of the day and with the output
  restricted to a finite set (F7); for F0, F2, F3a and F4 the bits
  `StartUpStabSingle`, `DeepByF1Rows` and `ReachRows` select the stability
  after a start-up by single instants, the deep decreases by the rows of the
  default formulation and the rows of the reach of the output from a
  start-up and towards a shut-down; every formulation is written on top of
  any formulation of `ThermalUnitBlock` and is kept up to date by its
  setters, or refuses the change

- `NuclearRules.h`: the codes of the formulations of the operating rules,
  their data, the labels of the states and the moves out of them, which
  `NuclearUnitExtDPSolver` and the flows on the label graph share

- the unit test checks every formulation of the operating rules against
  `NuclearUnitExtDPSolver` (the edge cases of the rules and random data),
  the chain of their continuous relaxations, the convex hull of the move of
  F2 and the setters after the generation (`NUBCfg-T-F*.txt`,
  `NUBCfg-T-TR-TM.txt`), that a `compute()` of `NuclearUnitExtDPSolver`
  that throws leaves the Solver usable, and its value after a change of
  `InitialPower` against a fresh unit

### Fixed

- the maximum power rows of the pt, SU, SD and SUSD formulations of
  `ThermalUnitBlock` (and so of `NuclearUnitBlock`) cap every instant of a
  run by the start-up limit plus the ramps up from the start-up and by the
  shut-down limit plus the ramps down to the shut-down, the first and the
  last instant included, where only the limits were used: the two arcs of a
  run that reaches the end of the horizon, which describe the same schedule,
  had different caps at its last instant, and the continuous relaxations of
  the pt, SU and SD formulations could be weaker than that of the DP one; a
  cap below the minimum power is no longer raised to it, so that a run that
  the limits and the ramps make impossible at an instant is excluded there

- `ACNetworkBlock` without `LineMinAngle` and `LineMaxAngle` no longer reads
  out of empty vectors when its rows are generated: the angle differences
  are then not bounded, there are no rows `AC_angle_bounds_limit` and
  `AC_elem_bounds`, and the strengthened relaxation, which needs the
  bounds, is not generated unless a Configuration asks for it, in which
  case it is refused; one of the two variables without the other is refused
  when read (`ACNetworkData::has_angle_bounds()` tells whether they are
  there)

- `ThermalUnitDPSolver` and `ThermalUnitExtDPSolver` compute again the
  ramps they take for a unit without `DeltaRampUp` or `DeltaRampDown` when
  `InitialPower` changes, so that the move from a new `InitialPower` above
  every `MaxPower` is free, as in the rows

- `ACNetworkBlock` reads `LineMinAngle` and `LineMaxAngle` as MATPOWER
  does: a value at most -360 (respectively at least 360), or NaN, leaves
  that side of the angle difference unbounded, and 0 for both leaves it
  unbounded altogether (`ACNetworkData::get_angle_difference_bounds()`);
  the bounds of +-360 of the MATPOWER instances made the rows
  `AC_angle_bounds_limit` force the imaginary part of the voltage product
  to 0 and `AC_elem_bounds` its real part above the product of the minimum
  voltages, and 0 / 0 made the angle difference 0. A line with an unbounded
  side, or with a range of a whole turn, now has neither group, and one
  with a range wider than 180 degrees has no `AC_angle_bounds_limit` (no
  such row is valid there); two finite bounds in the wrong order are
  refused when read

- the bounds `AC_elem_bounds` of `ACNetworkBlock` are the exact ranges of
  the real and imaginary parts of the voltage product over the bounds of
  the angle difference and of the voltages: the bound
  `sin( max - min )` of the imaginary part cut feasible points for
  symmetric bounds wider than +-45 degrees, and the lower bound of the real
  part was wrong for bounds beyond +-90; the rows `AC_angle_bounds_limit`
  of a bound beyond +-60 degrees are written as
  `cos( phi ) s - sin( phi ) c`, valid beyond +-90 too

- the strengthened relaxation of `ACNetworkBlock` covers the AC lines whose
  bounds on the angle difference are both within +-90 degrees, where its
  envelopes are valid, and leaves out the others, instead of using the
  envelopes of any range for every line

- the elementary check of the flows of `ACNetworkBlock` no longer reads the
  name of a line out of an empty vector when the data have no `LineName`

- the row of `NuclearUnitBlock` that forbids a modulation step at a start-up
  instant (`NoStartUpModulation_Nuclear`) is `m_t + v_t <= 1`, as documented,
  and no longer the void `m_t - v_t <= 1`, which let the MILP take a step at
  a start-up that the dynamic program does not admit

- `NuclearUnitBlock::is_feasible()` checks the part of `ThermalUnitBlock`
  with the tolerance it is given, and no longer with the default one of
  `Block`

- `NuclearUnitBlock::get_Solution()` saves the two auxiliary indicators of
  the deep decreases only if the unit has them

- `compute()` of `ThermalUnitDPSolver` and `ThermalUnitExtDPSolver` (and so
  of `NuclearUnitExtDPSolver`) releases the mutex of the Solver also when it
  throws, e.g., on a fixed Variable whose value the dynamic program refuses,
  and `NuclearUnitExtDPSolver::get_var_solution()` releases the lock of the
  Block likewise: the mutex stayed locked, and any other thread asking for
  the Solver waited forever

### Changed

- the instances of `data/` are downloaded and extracted by the targets
  `download_uc_<fmt>` and `extract_uc_<fmt>`, written as in every module
  that keeps its instances in the Package Registry, and the marker of the
  extraction carries the format in its name, so that a tree extracted
  before extracts once more

- the single-bus thermal instances (`tools/json2netCDF/instances_singlebus`
  and the `UC_singlebus` data) start from an initial state
  consistent with the demand and the reserve at instant 0: the units are
  committed in merit order while their minimum powers fit the lowest load net
  of the reserve of the first day, each holds an output with one ramp of room
  on both sides, and the others are off long enough to start at once; with the
  reserve, the secondary requirement is capped so that the reserve plus the
  change of the load in one instant fit within 0.7 times the one-instant ramp
  of the fleet; every thermal instance with the reserve is feasible over its
  horizon

- `json2nc4.jl` keeps the band breakpoints and the deep-decrease threshold of
  a `NuclearUnitBlock` at least a tenth of a ramp away from a whole number of
  ramps from the minimum and the maximum power, and the standalone units it
  writes enter the horizon on at their minimum power whatever the initial
  state of the fleet they are taken from, so that the standalone units of
  `1UC_Data/nuclear` differ from the previous ones only in those thresholds

- the ramp rows of the DP formulation of `ThermalUnitBlock` find the power
  of the previous instant of the run and the commitment of the run without
  scanning all the powers and all the runs for each row (the same rows)

- the minimum power, maximum power and initial perspective cut rows of the
  DP formulation of `ThermalUnitBlock` find the commitment and the cut
  variable of the run without scanning all the runs and all the cuts for
  each row (the same rows)

- `NuclearUnitExtDPSolver` takes the labels of its states and the moves out
  of them from `NuclearRules` (the same values and schedules as before); the
  bands and the deep decreases of `NuclearUnitBlock` are written by
  `build_band_rows()` and `build_deep_rows_f1()`, which the formulations of
  the rules share

## [0.10.0] - 2026-10-09

### Added

- `BatteryUnitBlock` has a kappa of its own for the converter
  (`ConverterKappa`, `set_converter_kappa()`, `get_converter_kappa()`,
  `get_converter_kappa_linearization()`), so that storage and converter can
  be sized apart; without it the converter follows `Kappa` as before

- the unit test checks the continuous relaxation of the SUSD formulation
  against that of the DP one (`LPRelaxBSCfg.txt`, `LPRelaxMILPCfg.txt`)
  and the two kappas of `BatteryUnitBlock`

- the unit test checks `update_rows()` on the rows of a class derived from
  `ThermalUnitBlock` and `set_init_updown_time()` after the generation in the
  seven formulations, and the harness of the power limits changes
  `InitUpDownTime` keeping and changing the initial state

- `ThermalUnitBlock` has the box of the active power, `0 <= p[ t ] <=` the
  operational maximum power, as the static group of `BoxConstraint`
  `ActivePowerBound_thermal`, kept up to date by `set_maximum_power()` and
  `set_availability()`: the rows of the commitment already imply it, but a
  Solver that only reads the boxes (say, a `BoxSolver` bounding the Objective,
  as `LagrangianDualSolver` does) needs it to see the power bounded

- the unit test checks the cyclic reservoirs of `HydroUnitBlock` after a
  change of the inflows and of the initial volume, the design variables of
  `BatteryUnitBlock` and `IntermittentUnitBlock` whose two bounds are 1, a
  battery of several modules with binary variables and reserves, and the
  refusal of a negative `NetworkCost`

- the unit test checks the costs and the rows at instant 0 of the reserves
  of `ThermalUnitBlock` against the DP Solvers and a :MILPSolver on random
  instances, the availability in the DP Solvers, the changes of the maximum
  power and of the availability, the reactive power of absorbing
  units, the bounds on the node injections after a scaling, the inertia rows
  of a unit at two nodes of a zone and the infinite storage bounds of the
  batteries

- the unit test compares, after a change of `MaxPower`, `Availability`,
  `InitialPower` or `InitUpDownTime`, the model of every formulation of a
  `ThermalUnitBlock` and of a `NuclearUnitBlock` with the one of the unit
  read afresh from the changed data, row group by row group, and the dynamic
  programming Solver with a :MILPSolver (`test_power_limits.cpp`, the
  formulations as `BlockConfig` files `TUBCfg-*.txt` and `NUBCfg-*.txt`; a
  reduced matrix by default, the whole one with `--power-limits --full`);
  it also checks the shut-down at instant 0 in every formulation and the
  setters of `ThermalUnitBlock` with an unordered `Subset`

- the page "The unit commitment model of UCBlock" (`UCBlockModel.h`) gives
  the notation, the conventions on instants, signs and scale factors, the
  composed model, the capacities and the investment, the Lagrangian dual
  with the meaning and sign of the dual values, the multistage use and one
  list of what is not modeled; the headers of every Block and DP Solver of
  the module state their rows as built, with one equation numbering each,
  and refer to it

- the unit test covers the network formulations on components, references
  and pure HVDC networks, the hydro delays, spillways, data checks and
  inertia, the storage balance, converter and design rows of the batteries,
  the design reserve rows of the intermittent units, the slack unit, the 3bin
  and DP maximum power rows and the initial power of `ThermalUnitBlock`, the
  end-of-horizon and band rules of `NuclearUnitBlock`, the zones and network
  constants of `UCBlock` and, when a :MILPSolver is built, the sign of the
  dual values of the linking rows

- data/gen_network_cases.py writes the instances of the edge cases of the
  network, small PyPSA dispatch networks whose DC lines have a susceptance
  (parallel lines in the same and in opposite directions, a triangle with a
  parallel side, DC components joined by a link, nodes reached only by
  links with losses), together with the optimum of PyPSA as their
  reference; they are in pypsa-data/ucblock of the data from 2026-10-03

- the module has a unit test of its own in `test/`, which needs nothing but
  the core SMS++ and builds all of its instances in memory: the three DP
  Solvers are compared with each other and with a brute force on small units
  (a horizon of one instant, minimum up and down times cut by the initial
  state or by the end of the horizon, ramps with the limits of the start-up
  and of the shut-down, the commitment fixed in some instants, the scale
  factor, the modulations of a nuclear unit), every UnitBlock and
  NetworkBlock and a `UCBlock` go through a netCDF round trip, and the
  setters of a `ThermalUnitBlock` and of a `UCBlock` are checked against a
  `FakeSolver`. The CI of the module builds it with `SMS++` alone and runs it

- the data archive 2026-09-26, which adds to `pypsa-data` the networks of
  pypsa2smspp that the batteries of `tests` read in the form of each module:
  one scenario of the modular family with the design in the units
  (`ucblock/smspp_mod_t48_s1_b2c_det_ucblock.nc`), the modular family as an
  MSSB (`mssb/smspp_mod_t48_s10_b2c_o2_mssb_ucblock.nc`) and a TSSB of the
  thermal family, whose scenarios are unit commitments
  (`tssb-thermal/smspp_tuc_u10_t24_s3_b1.nc4`)

- `ThermalUnitBlock::is_sol_feasible()` reads the schedule out of the
  `:UnitBlockSolution` it is given and checks it against the data of the unit,
  i.e. the operational bounds of the power, the reserves it takes room for, the
  ramps with the limits of the start-up and of the shut-down and the minimum up
  and down times with the state the unit comes from, so that the Variable of
  the unit are not needed, which `is_sol_feasible_physical()` says; a Variable
  of the schedule that is fixed (start-up, shut-down and reserves included)
  holds it to its value. Those are the constraints of the unit for an integral
  commitment, which is what every formulation of it encodes, hence a commitment
  that is not integral is not declared feasible; a unit that carries something
  the schedule does not answer for, i.e. a dimensioning variable, the reactive
  power, a reference schedule, a scale of its own or a fixed Variable of a
  formulation that the schedule only implies (the commitment differences, the
  pieces of the power, the cuts, the deviation from the reference schedule), is
  left to the check of the base class, which goes through the Variable and pays
  an allocation and two passes over them for each entry of a global pool that
  is revalidated; so is a `NuclearUnitBlock`, whose modulation the schedule
  does not show

- `tools/replicate_units`, which writes an instance with K times the thermal
  units of a given one, each copy carrying the data of the original, and
  multiplies by K what counts the replicated components, i.e., the
  generators, the storages, the demands and the budgets of the pollutants

- `ThermalUnitBlock::set_min_up_down_time()` sets the minimum up and down
  times before the Variable are generated, which is when they can still be
  set: unlike the other data they decide how many Variable there are, and the
  method throws once those exist

- the costs of the operating rules of a `NuclearUnitBlock`, i.e., those of a
  downward modulation step and of a deep decrease, can change:
  `set_down_modulation_costs()` and `set_deep_decrease_costs()` do it, each in
  its Range and Subset form and both in the methods factory, and a
  change of the corresponding coefficients of the Objective is folded into
  them, so that the physical representation follows and the DP Solver hears of
  it. They used to be refused, which stopped any Solver that writes on the
  Objective of the unit, such as the primal proximal heuristic of
  `LagrangianDualSolver` with its penalty term

- `ThermalUnitBlock` reads the optional variable `ShutDownCost`, the cost the
  unit pays at the instant it goes off, mapped over the time horizon as
  `StartUpCost` is and changed with `set_shutdown_costs()`. A unit that pays
  nothing to shut down keeps the vector empty and its Objective has no term
  for the shut-down variables, as before; the two DP Solvers charge the cost
  on the arc that closes an ON run

- the pollutant budget constraints of `UCBlock`: for each zone of each
  pollutant, the emission of the electrical generators at the nodes of the
  zone summed over the whole time horizon, i.e., the conversion factor
  `PollutantRho` (which may depend on time, and includes the duration of the
  time instant) times the active power, lies between `PollutantMinBudget` and
  `PollutantBudget`, the two matching a value when they are equal. The
  constraints are one vector with the zones of all the pollutants one after
  the other, as the budgets are in the file, they follow the scaling of the
  units, their duals are in the `UCBlockSolution` (bit 128), and the two
  budgets are changed with `set_pollutant_budget()` and
  `set_pollutant_min_budget()`; if `NumberPollutantZones` is not given every
  pollutant has one zone, and if `PollutantZones` is not given then all the
  nodes are in it. The constraints are grouped by pollutant, each with as many
  as its zones (`get_pollutant_constraints()[ p ][ z ]`), while the budgets
  and the duals keep the zones of all the pollutants one after the other, as
  the file does. They also take the levels of the storages into account
  (`PollutantStorageRho`, a factor on the level of each storage at each time,
  the storages of a unit being at the node of its first generator), which is
  how a limit that charges a storage for the change of its level is written;
  to this end `UnitBlock` has `get_number_storages()` and
  `get_storage_level()`, returning 0 and nullptr unless redefined, as
  `BatteryUnitBlock` (its charge), `HydroUnitBlock` (its reservoirs) and
  `HydroSystemUnitBlock` (the reservoirs of all its units) do

### Changed

- `BatteryUnitBlock::get_kappa_linearization()` leaves out the converter
  rows when the converter has a kappa of its own

- `ThermalUnitBlock::set_init_updown_time()` after the generation of the
  Variable takes a value that keeps whether the unit is on before the horizon
  and the first instant at which it may switch (any two values not smaller
  than `MinUpTime`, or any two not larger than `-MinDownTime`), the rows being
  written anew by `update_rows()`, and refuses any other one with the data
  kept, since it would make Variable appear or disappear

- `DCNetworkBlock` refuses a negative `NetworkCost`, when read and in
  `set_network_cost()`, since it makes the problem unbounded

- the documentation states the dual function of the version that relaxes
  the network rows, the nodal prices in the general PTDF formulation, the
  reactive multipliers of the Lagrangian dual, the derivations of the PTDF
  and CYCLE rows, of the battery binary rows and of the multistage cuts, and
  uses one notation in all the headers

- the costs of the reserves of `ThermalUnitBlock` are 0 when the data do not
  give them, `PrimaryRho` and `SecondaryRho` being no price, in the
  Objective and in the DP Solvers; a reserve with cost 0 is in the Objective
  only if the objective `Configuration` (bit 0 primary, bit 1 secondary)
  asks for it, and its cost cannot then become nonzero after the generation

- the reserve deliverability rows of `ThermalUnitBlock` are also written at
  instant 0 for a unit on before the horizon, from `InitialPower`, in every
  formulation and in the DP Solvers, and they follow `set_initial_power()`

- the DP Solvers of `ThermalUnitBlock` use the operational power bounds,
  `Availability` included, and without ramp data a ramp no move reaches

- `ThermalUnitBlock::set_maximum_power()` and `set_availability()` update every
  formulation once the constraints are generated: the rows are written in one
  place, `build_rows()`, which `update_rows()` runs again on the changed data
  and compares with the rows there are, changing coefficients and sides in one
  GroupModification, `NuclearUnitBlock` included; the new data are checked
  first (availability in [0, 1], `MaxPower` not below `MinPower`, a given
  `StartUpLimit` or `ShutDownLimit` within the operational bounds) and restored
  if the update fails; `MaxRampUpSteps` and `MaxRampDownSteps` in the data are
  refused; `update_rows()` refuses a change after which no row would be written
  to a static group that has rows (one written by `size_rows()` or
  `push_row()`, as the rows of a derived class are), or a row would be written
  to a dynamic group that is not generated

- the rows of `ThermalUnitBlock` whose number depends on the power limits
  are dynamic groups that the setters add to and remove from: the bound rows
  5 and 6 of the maximum power of the T formulation
  (`MaxPower5_Const_Thermal`, `MaxPower6_Const_Thermal`) and the ramp rows
  of several steps of the SUSD formulation (`RampUp_SUSD_Const_Thermal`,
  `RampDown_SUSD_Const_Thermal`); the rows of the T formulation whose
  Variable depended on the data have them all, with coefficient 0 where the
  term vanishes

- a `ThermalUnitBlock` on before the horizon shuts down at instant 0 only if
  `InitialPower` is at most `ShutDownLimit[ 0 ]`, whatever the ramps, in
  every formulation and in both DP Solvers (without `DeltaRampDown` the 3bin
  and T formulations and `ThermalUnitExtDPSolver` allowed it for any
  `InitialPower`): the interval ( 0 , 0 ) of the path formulations always
  exists and the row `ShutDownZero_Const_Thermal` fixes it (or `w_0` in
  3bin and T without `DeltaRampDown`) to 0 when it is not allowed

- `ThermalUnitBlock::set_initial_power()` updates every formulation once the
  constraints are generated, and no longer throws when `InitialPower` crosses
  `ShutDownLimit[ 0 ]` or the formulation is a path one with ramps

- a `StartUpLimit` or `ShutDownLimit` that the data do not give is the
  operational minimum power also after a change of the availability, so
  that an unavailable unit cannot produce

- the reactive power of `ThermalUnitBlock` is free in sign, and its bound
  rows, like those of `SlackUnitBlock`, always exist with the reactive power,
  an absent bound being 0

- the bounds on the node injections are given to the `NetworkBlock` again
  when the scale, the kappa or the maximum power of a unit changes, and they
  include kappa; the rows of the bounds follow them

- the KIRCHHOFF formulation of `DCNetworkBlock` fixes one voltage angle in
  each component of the lines with nonzero susceptance

- the Tikhonov coefficient of `DCNetworkBlock` is 0 also in the constructor
  and in `DCNetworkData::get_PTDF()`, and
  `DCNetworkBlock::change_active_demand_constraints()` is virtual

- `DCNetworkData::deserialize()` rejects a `ReferenceNode` that is not a node,
  and a network with hyperarcs and lines with nonzero susceptance

- a modulation of a `NuclearUnitBlock` with `PowerBands` crosses exactly one
  band boundary, the steps before the last keeping the output in the band of
  origin, in the rows and in `NuclearUnitExtDPSolver`

- `NuclearUnitBlock::set_initial_power()` throws when the new initial power
  is in another band than the one the rows at instant 0 were built for

- `ThermalUnitBlock` rejects a negative `FixedConsumption`, and its checks of
  `InitialPower` against the ramps at instant 0 throw once, with the name of
  the method

- `HydroUnitBlock` requires `LinearTerm` when some arc is a turbine, requires
  `StartArc` and `EndArc` together (and with more than one reservoir), checks
  their range and self-loops, and rejects an arc with `MinFlow` < 0 <
  `MaxFlow` or that is a turbine at some instant and a pump at another

- `HydroSystemUnitBlock::deserialize()` rejects a concave future cost of the
  water

- `BatteryUnitBlock` refuses a positive `MinPower` when intake and outtake are
  split, and `serialize()` writes `ConverterMaxPower` only if it was given

- `UCBlock::deserialize()` refuses `PrimaryZones`, `SecondaryZones` or
  `InertiaZones` absent with more than one zone of that kind, instead of
  reading out of bounds

- the `is_feasible()` of `ThermalUnitBlock`, `NuclearUnitBlock`,
  `HydroUnitBlock`, `HydroSystemUnitBlock`, `BatteryUnitBlock`,
  `IntermittentUnitBlock`, `SlackUnitBlock`, `DCNetworkBlock`,
  `DesignNetworkBlock`, `ECNetworkBlock` and `OTSNetworkBlock` accept the
  relative violation `Block::DefaultFeasTol` of the core when no
  Configuration gives a tolerance, instead of none

- the data archive is downloaded by version: `DATA_VERSION` in CMakeLists.txt
  names the version of the Package Registry to read, and the archive and the
  marker of its extraction carry it in their name, so that a tree holding
  an older extraction (the cache of the CI, or a clone extracted before)
  downloads and extracts again instead of running on the old data;
  data/upload-nc4 publishes the archive under that version

- the folders of the instances have one sub-folder per kind of problem, and
  the converters write into the one of the problem they translate, an
  instance of one kind having landed among those of another

- the fleet instances carry the rules that follow a start-up and a cycling
  copy of the load, so that the units they hold are exercised on the whole
  set of the operating rules and not on the part of it a flat load reaches

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- a `DCNetworkBlock` takes the KIRCHHOFF formulation where no Configuration
  says which one it wants, the PTDF one being both larger, since it carries a
  dense row per line, and the one whose flows a mixed network of AC lines and
  HVDC links used to get wrong; a `SimpleConfiguration< int >` of value 0 in
  the static-variables slot of the `BlockConfig` still asks for the PTDF
  formulation, and one of value 1 for the CYCLE one. A network whose lines all
  have a zero susceptance is a transport model under either formulation

- a `UCBlock` whose `NumberElectricalGenerators` is not the number of
  generators its units have is refused, instead of being read with the number
  the file states: everything indexed over the generators, from `GeneratorNode`
  to the emission rates, would be read over the wrong length, and the rows the
  `UCBlock` builds over them are sized with it, so that the model it gives
  depends on how far the two numbers are apart

- `ThermalUnitBlock` has the new virtual `update_objective_tail()`, with
  which a derived class makes the coefficients it appends to the Objective
  follow the scale factor: what includes its header has to be rebuilt

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

### Fixed

- the data check of a cyclic reservoir of `HydroUnitBlock` asks the total
  inflow over the horizon to lie between the smallest and the largest total net
  outflow that the flow bounds allow, a necessary condition, instead of the
  inflow of each instant to fit the arcs leaving the reservoir, which refused
  feasible instances; it is made also for a single reservoir without `StartArc`
  and `EndArc`

- the dynamic programming Solvers of `ThermalUnitBlock` no longer limit
  the move from `InitialPower`, when there is no `DeltaRampUp`/`DeltaRampDown`,
  by the largest `MaxPower`: a unit on before the horizon with `InitialPower`
  above `MaxPower` (e.g., one unavailable over the whole horizon) was declared
  infeasible, or given a power above its minimum at 0, while the formulations
  have no ramp row there

- the SUSD formulation of `ThermalUnitBlock` writes also the single-step
  ramp rows downwards on the runs by start-up and upwards on those by
  shut-down (`RampDown_Const_Thermal`, `RampUp_Const_Thermal`), as the SU
  and SD formulations do; without them its continuous relaxation could be
  weaker than that of the DP formulation

- the setters of the modulation ramps of `NuclearUnitBlock` change the rows
  through `update_rows()`, so that the coefficient of the modulation in the
  ramp rows with the direction of a modulation, the rows of TightRamp and the
  rows (32) of TightRules follow the new data; a change that adds or removes
  rows (32) is refused with the data kept, an unordered `Subset` keeps its
  values, and the `Range` form of `set_modulation_ramp_down()` issues
  `eSetModDM`

- `NuclearUnitExtDPSolver::load_fixings()` releases the read lock of the Block
  also when it throws

- `DCNetworkBlock` does not list "Kappa" among the expected netCDF names,
  since it does not read it

- the setters of `ThermalUnitBlock` with an unordered `Subset` sorted the
  indices but not the values, which went to the wrong instants (e.g.,
  `set_maximum_power()` with the values 60, 140 at the instants 4, 1 set 60
  at 1): the pairs are sorted together

- `NuclearUnitExtDPSolver` refuses a fixed band of the output, deep drop,
  deep low, modulation start or modulation end, whose fixings its labels
  cannot honor, instead of giving a schedule that may violate them, save the
  fixings at 0 of the first instant of a unit initially off;
  `NuclearUnitBlock` has the accessors of the band, modulation start and
  modulation end Variable

- `HydroUnitBlock::set_inflow()` and `set_initial_volume()` add the
  initial volume to the water balance of instant 0 only when it is
  nonnegative, as the generation does: a cyclic reservoir no longer gets
  the negative marker in its balance after a change of its inflows

- a design variable of `BatteryUnitBlock`, `IntermittentUnitBlock` or
  `DesignNetworkBlock` whose two bounds are 1 is fixed to 1 by its bound
  row; it was only bounded in [ 0 , 1 ], so the asset could be left out

- the binary rows and the reserve bounds of a `BatteryUnitBlock` whose battery
  design may exceed one module capped its power and reserves at those of one
  module: the binary rows are written with the largest design, and each reserve
  is bounded by the reserve of the modules built, a row with the design
  variable (`Primary_Design_Battery`, `Secondary_Design_Battery`), besides the
  bound of the largest design; an integer design may have a minimum above 1

- the cost of the secondary reserve of `ThermalUnitBlock` set after the
  generation landed on the wrong variable when the unit pays to shut down

- the inertia and reserve rows of `UCBlock` were updated with wrong
  positions for a unit with generators at several nodes of a zone

- `BatteryUnitBlock` with one instant and the cyclic level put the level
  twice in the row, and an infinite storage bound gave an infinite or NaN
  coefficient to the design variable (NaN bound with kappa 0)

- `ACNetworkBlock` gave the flows to the Solvers twice

- `ThermalUnitBlock` sized the rows fixing the commitment over the horizon
  instead of the instants it is fixed

- the PTDF and CYCLE formulations of `DCNetworkBlock` impose the balance of
  each component of the lines with nonzero susceptance, so that power does
  not move between components with no line carrying it, and `is_feasible()`
  checks these rows

- the PTDF formulation of `DCNetworkBlock` keeps the column of a
  `ReferenceNode` that is not the lowest node of its component

- `DCNetworkData::get_PTDF()` throws when the reduced susceptance matrix
  cannot be factorized, and without arguments it does not read the
  susceptances of a pure HVDC network out of bounds

- a change of the demand after generation updates the node balances of the
  CYCLE formulation of a pure HVDC network, does not write out of bounds in
  the KIRCHHOFF formulation of a mixed network, and updates the active
  balances of an `ACNetworkBlock` instead of crashing

- `OTSNetworkBlock` puts the cost of each flow in its Objective once

- `DCNetworkBlock::get_dual_prices()` gives a line with a design variable the
  sign of the dual of a bound

- a negative `UphillFlow` of `HydroUnitBlock` does not remove the withdrawal
  of the arc from the water balance of its start reservoir

- the flow-to-power row of a single-piece turbine of `HydroUnitBlock` is an
  equality only if a spillway with the same start and end reservoirs and the
  same delays can take its water at every instant, and the turbine has no
  positive `MinFlow` and no ramp rows

- `HydroUnitBlock::set_initial_volume()` throws if it would switch the cyclic
  closure of a reservoir on or off after the constraints are generated

- a `HydroSystemUnitBlock` without a PolyhedralFunctionBlock does not crash
  when its Variable are generated, and `get_polyhedral_function_block()`
  returns nullptr; `get_hydro_unit_block()` rejects the index
  `get_number_hydro_units()` in debug builds

- a change of the `InertiaPower` of a `HydroUnitBlock`, also inside a
  `HydroSystemUnitBlock`, changes the coefficients of the inertia rows of the
  `UCBlock`; adding an inertia term to a unit that had none when the rows
  were generated throws

- with only `ExtractingBatteryRho` given, the intake of a `BatteryUnitBlock`
  enters the storage balance with the sign of a charge

- `BatteryUnitBlock::set_initial_storage()` refuses, after the constraints are
  generated, a change of sign of `InitialStorage`, and `set_cost()` with an
  unsorted Subset gives each value to the index it comes with

- `ConverterMaxPower` bounds the converter of a `BatteryUnitBlock` also when
  the battery has no design variable, if the datum is given, and a converter
  design variable with no battery design variable is in its rows

- an `IntermittentUnitBlock` with a design variable has the reserve rows,
  with the capacity of the design variable in them

- `get_kappa_linearization()` of `BatteryUnitBlock` and
  `IntermittentUnitBlock` throws with a design variable instead of returning
  a wrong value

- `SlackUnitBlock` does not read `ActivePowerCost` out of range for the cost
  of the reactive power when the datum is absent, and does not throw in a
  `UCBlock` with inertia zones when the unit has no `MaxInertia`

- the 3bin maximum power rows of `ThermalUnitBlock` are counted and built
  right with a horizon of one instant or with the commitment fixed over the
  whole horizon, and with `MinUpTime` 1 a start-up at instant 0 is capped by
  `StartUpLimit` and the last power before a shut-down by the shut-down limit
  not above the maximum power, as in the other formulations and in the DP
  Solvers

- in the DP formulation of `ThermalUnitBlock` a run of one instant is capped
  by both the start-up and the shut-down limit

- the `InitialPower` of a `ThermalUnitBlock` on before the horizon is raised
  to `MinPower[0]`, also when the unit is unavailable at 0, before the
  numbers of ramp steps are derived from it; a value above `MaxPower[0]`
  gives a warning

- `ThermalUnitBlock::set_initial_power()` and `set_init_updown_time()`
  recompute the numbers of ramp steps from the initial power of the SUSD
  formulation

- a modulation of a `NuclearUnitBlock` cut by the end of the horizon cannot
  start upwards from the highest band or downwards from the lowest one, so
  that the MILP agrees with `NuclearUnitExtDPSolver`

- the operating-rule rows of a `NuclearUnitBlock` off before the horizon
  ignore `InitialPower`, and `set_solution()` does not end a modulation that
  covers a horizon shorter than `MaxModulationLength`, and takes a full-ramp
  step at the last instant that leaves the band of origin as the last step
  of a modulation

- with a single primary, secondary or inertia zone, a node whose zone index
  is not smaller than the number of zones is in no zone, as with several
  zones

- scaling a unit at a node in no reserve or inertia zone does not read and
  write past the rows of the zones

- the `ConstantTerm` of a `NetworkBlock` group of a `UCBlock` is not
  overwritten when `NetworkConstantTerms` is absent, and in the bus case the
  constant terms of the network are kept in the Objective of the `UCBlock`

- the overall balance of the PTDF formulation of DCNetworkBlock counts the
  losses of the HVDC lines, as the CYCLE one does: an HVDC line with
  an efficiency below 1 was forced to carry nothing, and a node reached
  only by such lines shed its whole demand

- the CYCLE formulation of DCNetworkBlock handles parallel DC lines: the
  spanning forest and the fundamental cycles are built on the lines rather
  than on the pairs of nodes, so that of two parallel lines one is in the
  forest and the other closes a cycle with it, while both used to be taken
  as tree edges and forced to carry the whole flow between their nodes (on
  a PyPSA network with two parallel lines the optimum was 5.6e-3 too high);
  also, an HVDC line with an efficiency other than 1 enters the flow of a
  tree edge with the efficiency at its end node, as in the nodal balance

- `NuclearUnitExtDPSolver` ignored a modulation (or a downward modulation,
  or a deep decrease) fixed to 1 at an instant in which the unit is off,
  since the fixings of those Variable were only checked on the moves of an
  on unit: a modulation at t needs the unit on at t and not starting up
  there (m_t <= u_t and m_t <= 1 - v_t), which is now what such a fixing
  imposes, and the DP no longer goes below the MILP when the tester fixes
  the modulations (TUDPS_FIXMOD)

- `BatteryUnitBlock::get_kappa_linearization()` reads the binary rows with
  the right sizes, sums the power rows under a reserve, and counts the intake
  and outtake bounds only on their upper side

- the cost setters of a `ThermalUnitBlock` (`set_linear_term()`,
  `set_quad_term()`, `set_const_term()`, `set_startup_costs()`,
  `set_shutdown_costs()`, `set_reactive_linear_term()` and those of the
  spinning reserves) write in the Objective the scale factor times the cost,
  as the Objective is generated. They wrote the cost of one copy, so that on
  a unit with a `Scale` other than 1 the Objective no longer agreed with the
  data after a change

- the minimum up (down) time of a `ThermalUnitBlock` is bounded by the
  horizon plus the number of instants the unit has been on (off) before it,
  both when it is read and when `set_min_up_down_time()` sets it, so that
  the unit stays as it is for as long as the minimum time says. The bound was
  the horizon plus one, which freed a unit that had been on (off) for more
  than one instant before the minimum time was over: a unit off since 3
  instants with a `MinDownTime` of 4 could start at the first instant of a
  horizon of one

- `ThermalUnitDPSolver.h` can be included where FastFlow is not seen: the
  constructor, which needs the destructor of the `ff::ParallelFor` should it
  throw, is defined in the `.cpp` as the destructor is

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- the scaling of a unit rewrites the rows of that unit also while a
  `LagrangianDualSolver` is attached: the index of the unit was asked to its
  father, which is then the `LagBFunction` holding the unit alone, so that
  every unit was taken for unit 0; the index is now looked up among the
  units of the `UCBlock`

- the reserve band of a unit that is on for a single period, i.e., starts up
  at t and shuts down at t+1, is capped by the smaller of the start-up and the
  shut-down limit in both `ThermalUnitDPSolver` and `ThermalUnitExtDPSolver`,
  in their solutions as in their values, and so is that of a unit that is on
  before the horizon and shuts down at the first instant; the cap was that of
  the start-up alone, and the two DPs could give a value below the optimum

- `ThermalUnitExtDPSolver` builds the solution with the shut-down cap only
  where it bites, i.e., where it is below the maximum power, and restricts the
  value function to the domain of the shut-down before adding its term at the
  first instant, which was read at a point where the band was still uncapped

- the formulations of `ThermalUnitBlock` with time-varying minimum and maximum
  power are exact: in the T formulation the ramp rows weigh the state of the
  previous instant with the minimum power of the right instant, every
  coefficient of the max-power rows is measured against the maximum power of
  the row, and a single on period is capped by the smaller of the start-up and
  shut-down limits (the sign test was inverted and took the larger); in the
  pt, DP, SU, SD and SUSD formulations every shut-down cap is the limit of the
  instant the unit shuts down at, and the start-up limit of the SD ramp-up no
  longer reads one instant past the horizon

- in the design MIP of a unit with reactive power, the off-bounds of the
  reactive power are multiplied by the design variable, so that a unit that is
  not built produces no reactive power, as in the DPs

- `ThermalUnitDPSolver` with time-varying minimum and maximum power no longer
  writes out of bounds in the energy sweep when the lower end of the domain
  lies more than a ramp-down below the previous one, and a domain that is
  empty or below the shut-down cap makes the run infeasible instead

- the dominance test of `ThermalUnitExtDPSolver` between two value functions
  is linear in their pieces, and the profiling switches of the solver are read
  from the environment only when `TUEDPS_PROFILE` is set

- `ThermalUnitExtDPSolver` (and so `NuclearUnitExtDPSolver`) reads the label
  of the initial state as an off-state through `shut_label()` when the unit
  shuts down at the first instant or restarts from an initial off period,
  since the on-code of that label and the off-code a shut-down produces are
  different in general; the label was used as it was, i.e., as an on-code,
  and a unit whose initial state has no off-code does not take those arcs

- the PTDF matrix of a `DCNetworkBlock` is no longer perturbed by a Tikhonov
  term on the diagonal (`f_tikhonov_coeff` is 0 by default), the reduced
  Laplacian being nonsingular once a reference node per connected component is
  removed. The perturbation broke the conservation of the flows: the PTDF rows,
  the nodal balance of the nodes an HVDC line touches and the overall balance
  are then no longer consistent, and together they pin a relation among the
  injections that has nothing to do with the network, so that on a network of
  both kinds of lines the PTDF formulation gave an optimum far from that of the
  KIRCHHOFF one, and could even declare the problem infeasible. The two now
  agree, on a network of 8 countries as on a triangle

- two lines that join the same two nodes sum their susceptances in the PTDF
  matrix, instead of the second one overwriting the entry of the first

- `ThermalUnitBlock` accepts a `MinUpTime` (`MinDownTime`) of one instant more
  than the horizon, which says that the unit, being on (off) before it, never
  switches within it: the commitment is then fixed at every instant and there
  is no start-up or shut-down variable. The bound was the horizon itself, so
  that the last instant was free however long the minimum time

- the step that fetches the data archive of this module says what went wrong
  when it goes wrong: the download is checked, an archive that did not arrive
  is removed instead of being left on disk for the build to take for the real
  one, and the message names the URL. A server that answers with an error page
  used to leave a file of a few bytes there, which made the next build fail
  while extracting it, with the message of `tar` and no mention of the
  download

- scaling a `ThermalUnitBlock` after its Objective has been generated rewrites
  every term that carries the scale factor, i.e., also those of the shut-down,
  of the primary and secondary reserves, of the perspective cuts and of the
  reactive power, and in a `NuclearUnitBlock` those of the downward modulation
  steps and of the deep decreases, which used to keep the old factor; the model
  of a unit scaled after being built is now the same as that of a unit scaled
  before

- a `BatteryUnitBlock` reads the `ReferenceSchedule` its file declares: the
  variable was among the ones it expects and all the machinery was there, the
  deviation variables, the rows and the term of the Objective, but the datum
  was never deserialized, so that the schedule of a battery was thrown away in
  silence while a thermal or a hydro unit followed the one it is given. The AC
  instances carry one on each of their 15 batteries, as they do on their 153
  thermal units, hence the two kinds of unit now behave the same way

- a `BatteryUnitBlock` that follows a reference schedule refuses to have its
  kappa changed, i.e. to be resized: whether the profile a resized battery is
  asked to follow is the same one in absolute terms or one resized with it is
  not settled, and leaving the profile where it is would answer it in silence.
  `get_reference_schedule()` gives the schedule, empty where there is none

- the reactive node injection constraints of a `UCBlock` carry the reactive
  power of the generators and nothing else, as their documentation says: they
  used to carry the fixed consumption too, which is an ACTIVE power, on the
  commitment variables. Whether a unit that is off also absorbs reactive
  power is not settled; were it to, it would call for a datum of its own

- the fixed consumption of a unit that is off raises the right-hand side of
  the node injection constraints at a single node, as it already did with
  more than one node and as the setter of the demand already recomputed it:
  the generation lowered it instead, so that a unit consuming while off made
  the others generate less, and merely writing the demand back into the
  `UCBlock` changed the model. The value of the plan4res instances that carry
  a fixed consumption moves by 8 to 12%. The formula in the documentation,
  which wrote the consumption as a positive term of the injection, follows

- an `ACNetworkBlock` and an `OTSNetworkBlock` refuse a kappa on one of
  their lines, `DCNetworkBlock::set_kappa()` being virtual now and their
  override throwing: they build their own rows and not the power flow limit
  ones the kappa is written into, so that sizing a line of theirs used to
  reach rows that are not there. What supporting it would take is written
  where they refuse it. The method that writes the kappa into the rows
  refuses as well when there are none, so that any other derived class is
  told rather than left to write where there is nothing

- the dynamic programming Solvers refuse a unit that has a reference
  schedule, of which they have no term: they used to answer for a unit that
  pays nothing to depart from its schedule, i.e., a value that is not the one
  of the Objective, and silently. `ThermalUnitBlock::get_reference_schedule()`
  gives the schedule, empty where there is none

- the deviation from the reference schedule of a `ThermalUnitBlock` or of a
  `BatteryUnitBlock` is weighed with the scale factor, as every other term of
  their Objective: the schedule is that of one unit, from which each of the
  copies deviates on its own, so that the fleet pays the scale factor times
  what one copy pays; it used to be weighed with 1, so that the Objective of
  a scaled unit was neither the cost of one copy nor that of all of them. The
  guard that refuses a change of those coefficients now asks them to be the
  scale factor. A `HydroUnitBlock` has no scale factor, hence its own
  deviation is unchanged, and the reference schedule of a `BatteryUnitBlock`
  whose kappa changes is left as it is [see `set_kappa()`]

- a change of the coefficients of the Objective of a `ThermalUnitBlock`, as a
  dualizing Solver makes, is divided by the scale factor before being stored in
  the costs of the unit, which are those of one copy; they used to be stored as
  they are, i.e., multiplied by the scale factor, and a scaling made with
  `eModBlck` went the same way through the start-up and fixed costs

- the dynamic programming Solvers of the `ThermalUnitBlock` and of the
  `NuclearUnitBlock` report the value of all the copies of the unit, i.e.,
  the scale factor times that of the one copy they solve for, which is the
  value of the Objective; the Solution they produce is still that of one
  copy, as the Variable of the unit are

- the reactive node injection constraints of a `UCBlock` follow the scale
  factor of a unit as the active ones do, and a scaled unit with a fixed
  consumption at a single node gets the fixed consumption, and no longer the
  scale factor, as the coefficient of its commitment

- scaling a unit of a `UCBlock` with a single primary, secondary or inertia
  zone no longer reads the zone of its generators out of an empty vector,
  which crashed

- `IntermittentUnitBlock::scale()` tells the `UCBlock` even when no Solver
  is listening, as the thermal unit and the battery do, so that the rows
  carrying the scale factor are rewritten anyway

- the ramps of a `ThermalUnitBlock` that change over time are read as the
  documentation says in every formulation and Solver: `DeltaRampUp[ t ]`
  (`DeltaRampDown[ t ]`) bounds the increase (decrease) of the power from
  t - 1 to t, from `InitialPower` if t = 0. The 3bin and pt formulations,
  the reserve deliverability rows and the three dynamic programming Solvers
  used `DeltaRampUp[ t - 1 ]` for that step, so that `DeltaRampUp[ 0 ]` bound
  both the first two steps and the formulations disagreed with each other;
  the T formulation mixed the two indices in its ramp-down rows, and the
  formulations whose rows span several steps (SUSD, the bounds of T, the
  maximum powers of the interval formulations) multiplied one ramp by the
  number of steps instead of summing the ramps of the steps. With constant
  ramps nothing changes

- `IntermittentUnitBlock::check_data_consistency()` refused a
  `MinCapacityDesign` above 1 when `MaxCapacityDesign` is negative, as if the
  design were binary, while that design is an integer between 0 and
  `|MaxCapacityDesign|` (e.g., the number of modules of a modular asset): it
  now asks only that `MinCapacityDesign` be at most `|MaxCapacityDesign|`

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

- the band of a `NuclearUnitBlock` could move against the direction of the
  modulation that moves it, the rows asking only that it change by one: where
  the landing power is a breakpoint, and both bands hold it, a modulation
  upwards could therefore leave the unit a band lower, which costs nothing
  while a downward one is paid for. The move now follows the direction, as it
  does in the dynamic program, the two having disagreed by up to `1.2e-4` on
  the instances of the battery. `has_modulation_direction()` is true whenever
  the output is banded, the direction Variable being what says which way the
  band goes

- at the last instant of the horizon a `NuclearUnitBlock` with the bands could
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

- the cuts that `ThermalUnitBlock` and `NuclearUnitBlock` separate carried a
  constant term of 1 in their `LinearFunction`, `eNoMod` having been passed to
  its constructor, where it is the constant, rather than to `set_function()`.
  The `MILPSolver` reads the coefficients and not the constant, hence the rows
  it got were right, but every check done on the `Block`, `is_feasible()`
  comprised, saw each separated cut violated by 1 at any solution

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

### Removed

- test/ moved to ThermalUnitBlock_Solver in test repository

- useless test_package

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

[Unreleased]: https://gitlab.com/smspp/ucblock/-/compare/0.10.0...develop
[0.10.0]: https://gitlab.com/smspp/ucblock/-/compare/0.9.0...0.10.0
[0.9.0]: https://gitlab.com/smspp/ucblock/-/compare/0.8.0...0.9.0
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
