# json2netCDF

A small Julia tool that produces SMS++ `UCBlock` `netCDF` instances of the
**single-bus thermal Unit Commitment** problem. It takes randomly generated
realistic single-bus thermal UC instances, possibly with primary and secondary
reserve and reactive power constraints, in JSON format and converts them in the
SMS++ `netCDF` file consumed by `UCBlock`.

It is the JSON counterpart of the `csv2netCDF` tool (which targets the energy
community model): same layout, same `gen-all-nc4` driver, and the same idea of
committing only the lightweight source instances (here the JSON) while the
generated `.nc4` go to the standard `data/nc4/` location and are not committed.

## Layout

```
json2netCDF/
├── json2nc4.jl                 # the JSON -> SMS++ netCDF converter
├── gen-all-nc4                 # bash driver: (re)convert every JSON -> .nc4
├── gen-nuclear             # bash driver: the standalone nuclear units of
│                               #   the tests of the nuclear operating rules
├── gen-nuclear-scaling         # bash driver: the same units at the day, week
│                               #   and month horizons
├── gen-nuclear-uc              # bash driver: the single-bus UCBlock whose
│                               #   units are all nuclear
├── Project.toml / Manifest.toml # pinned env for the converter (JSON+NCDatasets)
└── instances_singlebus/        # the JSON source instances (sb_*.json)
```

The committed `Project.toml` + `Manifest.toml` pin only the converter's
dependencies (`JSON`, `NCDatasets`), so `julia --project=.` instantiates a
lightweight environment with no MILP solver.

The JSON under `instances_singlebus/` are the committed source of the suite;
the converter writes the `.nc4` into `../../data/nc4/UC_singlebus/`, the
standard location the tests and tools expect (and where the suite is
distributed from, not committed). The number of periods is therefore

| horizon | days | periods (= days × 96) |
|---------|------|-----------------------|
| `day`   | 1    | 96                    |
| `week`  | 7    | 672                   |
| `month` | 31   | 2976                  |

(note `month` is 31 days, i.e. 2976 periods, not 30.)

## The converter

`json2nc4.jl` reads a single-bus thermal UC JSON instance (schema
`v4.0-singlebus-reserve`) and writes an `UCBlock` `netCDF` with:

- one `ThermalUnitBlock` per generator (all on node 0);
- `ActivePowerDemand` on the single node;
- optional primary (FCR) and secondary (aFRR) reserve: a single reserve zone,
  `PrimaryDemand` / `SecondaryDemand` at `UCBlock` level and `PrimaryRho` /
  `SecondaryRho` at `ThermalUnitBlock` level.

All quantities are assumed already discretised to the instance time-step
(min up/down times in periods, ramps in MW/period, costs per period) and are
copied verbatim.

### Usage

```bash
# convert one instance to an explicit output
julia --project=. json2nc4.jl in.json out.nc4
# (re)convert the whole suite: every instances_singlebus/*.json -> data/nc4/UC_singlebus
./gen-all-nc4
```

The conversion is deterministic: `./gen-all-nc4` regenerates the suite `.nc4`
byte for byte under `data/nc4/UC_singlebus/`.

### Standalone units

With `--thermal-single <dir>` or `--nuclear-single <dir>` the converter writes,
instead of one `UCBlock`, one standalone `ThermalUnitBlock` or
`NuclearUnitBlock` file per generator (the first `--n-units`), which the
`ThermalUnitBlock_Solver` tester reads directly. A nuclear unit gets a
time-varying linear cost that follows the demand, so that it has a reason to
change its output, and the modulation data (`--mod-time`, `--mod-frac`). The
operating rules of `NuclearUnitBlock` are added by further options (see
`julia json2nc4.jl` with no argument for the list): the longest modulation,
the length of a day and the daily limits on the modulations and the
start-ups, the deep decreases and their daily limit, the costs of the
downward modulation steps and of the deep decreases, spinning reserves, a
commitment-gated reactive box and the number of periods. Each of them is only
written when it differs from its default, so that the default output does not
change. `./gen-nuclear` builds with them the instances of the test battery
of the nuclear operating rules under `data/nc4/1UC_Data/nuclear/48/`, and
`./gen-nuclear-scaling` the same units with the full set of the rules at the
horizons of one day, one week and one month, under
`data/nc4/1UC_Data/nuclear/{96,672,2976}/`.

### Nuclear UCBlock

With `--nuclear` the converter writes one `UCBlock` whose units are all
`NuclearUnitBlock`, with the same operating-rule options as above;
`./gen-nuclear-uc` builds the single-bus instances of the 50-unit fleets at
the three horizons under `data/nc4/UC_singlebus-nuclear/`.

## Authors

- **Antonio Frangioni**
- Dipartimento di Informatica, Universita' di Pisa
