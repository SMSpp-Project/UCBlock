# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- The way a Block is sized has a name in the methods factory, so that a
  consumer can write it and read the sensitivity back without knowing which
  class it is holding. `ThermalUnitBlock`, `NuclearUnitBlock`,
  `IntermittentUnitBlock` and `BatteryUnitBlock` register `scale()` as
  `<class>::replicate`, which stands for k identical copies of the unit;
  `IntermittentUnitBlock`, `BatteryUnitBlock` and `DCNetworkBlock` register
  `set_kappa()` as `<class>::resize`, which stands for one copy of k times the
  size. The thermal and the nuclear units deliberately do not offer the
  second road, their commitment being binary. The old names stay registered, a registered name travelling
  inside the instances that name it.

- The sensitivity that mirrors each of those names. `<class>::get_resize_
  linearization` gives the derivative of the value with respect to kappa, and
  is answered by the Block that owns the rows kappa writes into;
  `UCBlock::get_replicate_linearization` gives the derivative with respect to
  the number of copies of one of its units, and is answered by the container,
  because that factor lives in the rows the container builds on top of the
  unit and not in the rows of the unit. Both come in the range and the subset
  form, both refuse a span shorter than what they are asked for, and both
  refuse an index that does not exist. The derivative with respect to the
  copies sums, row by row, what the factor multiplies there: the power of the
  unit in the active and in the reactive node injection, with the fixed
  consumption of a unit that is off subtracted, its reserves, its inertia,
  its emissions in the pollutant budget with the levels of its storages, and
  the cost of one copy.

- `DCNetworkBlock::get_resize_linearization( line )` and
  `UCBlock::get_replicate_linearization( unit )`, the two protected methods
  those sensitivities are read from.

- the data archive 2026-09-26, which adds to `pypsa-data` the networks of
  pypsa2smspp that the batteries of `tests` read in the form of each module:
  one scenario of the modular family with the design in the units
  (`ucblock/smspp_mod_t48_s1_b2c_det_ucblock.nc`), the modular family as an
  MSSB (`mssb/smspp_mod_t48_s10_b2c_o2_mssb_ucblock.nc`) and a TSSB of the
  thermal family, whose scenarios are unit commitments
  (`tssb-thermal/smspp_tuc_u10_t24_s3_b1.nc4`)

- `ThermalUnitBlock::is_sol_feasible()` reads the schedule out of the
  `:UnitBlockSolution` it is given and checks it against the data of the unit,
  i.e. the operational bounds of the power, the reserves it takes room for,
  the ramps with the limits of the start-up and of the shut-down and the
  minimum up and down times with the state the unit comes from, so that the
  Variable of the unit are neither needed nor touched, which
  `is_sol_feasible_physical()` says. Those are the constraints of the unit for
  an integral commitment, which is what every formulation of it encodes,
  hence a commitment that is not integral is not declared feasible; a unit
  that carries something the schedule does not answer for, i.e. a
  dimensioning variable, the reactive power, a reference schedule, a scale of
  its own or a Variable that is fixed, is left to the check of the base class,
  which goes through the Variable and pays an allocation and two passes over
  them for each entry of a global pool that is revalidated

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

- `BatteryUnitBlock::get_kappa_linearization()` reads the two rows that keep
  a battery from charging and discharging at once with the sizes they carry,
  the minimum power on the charging row and the maximum one on the
  discharging row: it had them swapped, which gave a wrong sensitivity for a
  battery whose two powers differ, while a symmetric one hid it

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

[Unreleased]: https://gitlab.com/smspp/ucblock/-/compare/0.9.0...develop
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
