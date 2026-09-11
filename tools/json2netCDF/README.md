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

## Authors

- **Antonio Frangioni**
- Dipartimento di Informatica, Universita' di Pisa
