# csv2netCDF

This is a tool to convert UC data in the Energy Communities setting from
csv to NetCDF.

## Getting started

The main script file is `csv2nc4.jl` that optionally takes in input the
following parameters:

```sh
julia csv2nc4.jl [yml]
```

where `yml` can be one of the followings:

- `energy_community_model_CO`     (i.e., Cooperative case; the default if none is given)
- `energy_community_model_NC`     (i.e., No Cooperative case)

The same names with the `_sto` suffix enable scenario sampling (the
deterministic UCBlock is still produced; an additional `TSSB_*.nc4` is
written that references it via the `filename` attribute).

By default the thermal generator is included (the output file gets the `_TUB`
suffix and SMS++ models the generator as a `ThermalUnitBlock`). Pass
`--no-thermal` to opt out — same naming convention as the EC.jl@stochastic
flag of the same name:

```sh
julia csv2nc4.jl [yml] --no-thermal
```

The "No Asset" (NA) family of instances — previously generated from dedicated
`*_NA*.yml` files, now removed because they were identical to CO with all
installable assets disabled — is produced by passing `--no-asset` to the CO
YAML. The output is written as `*_NA_Test*.nc4`:

```sh
julia csv2nc4.jl energy_community_model_CO --no-asset
```

`--no-asset` implies no thermal (so it cannot be combined with `--no-thermal`).

For test purposes, the script can take an additional flag to enforce
the generation of the physical ECNetworkBlock(s), which adds the `_NB` suffix:

```sh
julia csv2nc4.jl [yml] --with-network-blocks
```

Finally, for `PV` / `wind` (`IntermittentUnitBlock`) and `batt` / `conv`
(`BatteryUnitBlock`) installable assets, the fleet of `N = max_capacity /
nom_capacity` identical modules can be encoded in three LP-equivalent ways
selected by `--design-mode`:

```sh
julia csv2nc4.jl [yml] --design-mode=fleet
julia csv2nc4.jl [yml] --design-mode=scale
julia csv2nc4.jl [yml] --design-mode=design   # default
```

- `fleet`: a single block sized by `max_capacity`; the design variable is
  continuous in `[0, 1]`, no `Scale` / `MaxCapacityDesign` emitted.
- `scale`: a single block sized by `nom_capacity` with `Scale = N`; the
  `f_scale` factor multiplies the cost and power coefficients so the block
  behaves as a fleet of `N` identical modules.
- `design`: a single block sized by `nom_capacity` with `MaxCapacityDesign
  = ±N` (`BatteryMaxCapacityDesign` / `ConverterMaxCapacityDesign` for
  batteries). The sign is negative ⇔ integer install ∈ `{0,…,N}` (chosen
  when the YAML asset has `modularity: false`, and always in the stochastic
  flow to match `EnergyCommunity.jl@stochastic`); otherwise continuous
  install in `[0, N]`.

The three modes are mathematically equivalent at LP-relaxation level; only
`design` with a negative `MaxCapacityDesign` actually enforces an integer
install at MILP level. `Scale ≠ 1` and `|MaxCapacityDesign| > 1` (or
`|Battery/ConverterMaxCapacityDesign| > 1` for batteries) are *mutually
exclusive* — `check_data_consistency()` on the C++ side rejects
configurations that activate both. Thermal assets are unaffected: a
`ThermalUnitBlock` has a binary design variable only, so a granular integer
count `{0,…,N}` is achieved by replicating the block `N` times rather than
via `Scale = N` (which would force a synchronous all-or-nothing fleet).

The convenience wrapper `gen-all-nc4` runs the deterministic and stochastic
batches over CO and NC (with `_TUB` default and `--no-thermal` variants) plus
the NA case via CO + `--no-asset`, with and without `--with-network-blocks`.

## YAML layout (single format)

- EC-wide profiles (`time_res`, `energy_weight`, `reward_price`,
  `peak_categories`) live under `general.profile`.
- Pricing fields (`buy_price`, `sell_price`, `consumption_price`,
  `peak_tariff`, `peak_weight`) live under per-tariff blocks selected by
  `users.<u>.tariff_name` (e.g. `commercial`, `non_commercial`). Per-user
  contributions are summed across users, so EC-wide prices scale with the
  number of users in `general.user_set`.

## Output

The script writes one or two files into `../../data/nc4/EC_Data/`:

- `EC_<MODE>_Test[_TUB][_NB].nc4` — the deterministic UCBlock instance, always
  produced.
- `TSSB_EC_<MODE>_Test[_TUB][_NB].nc4` — the TwoStageStochasticBlock instance,
  produced only when `scen_s_sample * scen_eps_sample > 1`. The TSSB does NOT
  embed the deterministic UCBlock inline; instead its inner `Block` group
  references the `EC_*.nc4` companion file via the `filename` attribute.

Stochastic scenario sampling uses `pem_extraction`, `scenario_definition` and
`Scen_eps_sampler`. Profiles whose `std` field is missing fall back to a
multiplicative noise `sigma * |mean|`.

## Author

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa
