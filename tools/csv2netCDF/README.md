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

The same names with the `_two_stage` suffix enable scenario sampling with a
single short-period scenario per long-period one, and the `_three_stage`
suffix enables several short-period scenarios per long-period one (the
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

For a stochastic YAML, `--multistage` builds a `MultiStageStochasticBlock`
that aggregates one `TwoStageStochasticBlock` per long-period (`scen_s`)
scenario and ties the first-stage variables across them, writing an
`MSSB_EC_*.nc4` instead of `TSSB_EC_*.nc4`:

```sh
julia csv2nc4.jl [yml]_two_stage  --multistage    # MSSB_EC_*.nc4
julia csv2nc4.jl [yml]_three_stage --multistage   # MSSB_EC_*_3S.nc4
```

With a `_two_stage` YAML (a single short-period scenario) the extensive form
coincides with the flattened `TwoStageStochasticBlock`. With a `_three_stage`
YAML (several short-period scenarios) the inner `TwoStageStochasticBlock` also
ties the day-ahead declared-dispatch bid across the short-period scenarios, so
the output carries an extra `_3S` suffix. The flag has no effect on a
deterministic YAML.

Finally, for `PV` / `wind` (`IntermittentUnitBlock`), `batt` / `conv`
(`BatteryUnitBlock`) and `generator` (`ThermalUnitBlock`) installable
assets, the fleet of `N = max_capacity / nom_capacity` identical modules
can be encoded in three LP-equivalent ways selected by `--design-mode`:

```sh
julia csv2nc4.jl [yml] --design-mode=fleet
julia csv2nc4.jl [yml] --design-mode=scale
julia csv2nc4.jl [yml] --design-mode=design   # default
```

- `fleet`: a single block sized by `max_capacity`; the design variable is
  continuous in `[0, 1]`, no `Scale` / `MaxCapacityDesign` emitted. For
  thermal, this produces a single `ThermalUnitBlock` with binary design
  ⇒ synchronous fleet (install ∈ {0, N}) — mathematically equivalent to
  `scale` mode for thermal.
- `scale`: a single block sized by `nom_capacity` with `Scale = N`; the
  `f_scale` factor multiplies the cost and power coefficients so the block
  behaves as a fleet of `N` identical modules. For thermal the commitment
  `u_t` is binary and shared by the `N` modules ⇒ synchronous on/off
  fleet (install ∈ {0, N}).
- `design`: a single block sized by `nom_capacity` with `MaxCapacityDesign
  = ±N` (`BatteryMaxCapacityDesign` / `ConverterMaxCapacityDesign` for
  batteries) for PV / wind / batt / conv. For thermal the block is
  instead **replicated `N` times** (each replica with its own binary
  design and independent `u_t`), since `ThermalUnitBlock` has no
  `MaxCapacityDesign` and the granular integer install ∈ {0,…,N} with
  *independent* commitments is achievable only by replication. The sign
  of `MaxCapacityDesign` (PV / wind / batt / conv only) is set by the
  per-asset YAML field `modularity`, using the same convention as
  `EnergyCommunity.jl`:

  | YAML `modularity`            | `MaxCapacityDesign` |
  | ---------------------------- | ------------------- |
  | `true`                       | `−N` (integer)      |
  | `false` (default if missing) | `+N` (continuous)   |

  The convention is applied uniformly on both the deterministic and the
  stochastic flow. EC.jl's stochastic branch currently keeps the
  first-stage `n_us` integer regardless of `modularity` (via
  `StochasticPrograms`), so for now an integer-flagged YAML matches
  EC.jl@stochastic exactly while a continuous-flagged YAML diverges
  there; honouring the field on the stochastic side too makes the
  encoding forward-compatible with a future EC.jl@stochastic that
  exposes the same toggle. The `modularity` field is consulted only in
  `design` mode: `fleet` and `scale` never emit a `MaxCapacityDesign`
  and always keep the design continuous, so for them the YAML field has
  no effect.

The three modes are mathematically equivalent at LP-relaxation level. At
MILP level: for PV / wind / batt / conv, only `design` with a negative
`MaxCapacityDesign` enforces an integer install; for thermal, only
`design` (replication) gives a strictly larger feasible set than `scale`
/ `fleet` thanks to independent commitments per replica. `Scale ≠ 1` and
`|MaxCapacityDesign| > 1` (or `|Battery/ConverterMaxCapacityDesign| > 1`
for batteries) are *mutually exclusive* — `check_data_consistency()` on
the C++ side rejects configurations that activate both.

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

## Resizing the instance

To change the instance size you only edit the YAML — no manual edits to the
market CSV are needed:

- **Fewer users**: shrink `general.user_set` to the users you want.
- **Fewer time steps**: set `general.final_step` to the desired value. On read,
  `read_input` recalibrates the market dataset so the annual energy balance
  `final_step * time_res * energy_weight = 8760` (hours in a year) is preserved:
  `time_res` (the physical step length, e.g. 15 min) is kept as-is and the
  annual repetition weight `energy_weight` is recomputed as
  `8760 / (final_step * time_res)`. The CSV is rewritten in place only when the
  value actually changes (the calibrated default `final_step = 96`,
  `time_res = 0.25`, `energy_weight = 365` is a no-op). Make sure
  `general.optional_datasets` points at the matching market CSV (the `*_sto`
  variants for the stochastic flow).

## Output

The script writes one or two files into `../../data/nc4/EC_Data/`:

- `EC_<MODE>_Test[_TUB][_NB].nc4` — the deterministic UCBlock instance, always
  produced.
- `TSSB_EC_<MODE>_Test[_TUB][_NB].nc4` — the TwoStageStochasticBlock instance,
  produced only when `scen_s_sample * scen_eps_sample > 1`. The TSSB does NOT
  embed the deterministic UCBlock inline; instead its inner `Block` group
  references the `EC_*.nc4` companion file via the `filename` attribute.
- `MSSB_EC_<MODE>_Test[_TUB][_NB][_3S].nc4` — the MultiStageStochasticBlock
  instance, produced with `--multistage` from a stochastic YAML. It aggregates
  the `TwoStageStochasticBlock` as its inner sub-Block(s). The `_3S` suffix
  marks a three-stage instance (from a `_three_stage` YAML), where the inner
  `TwoStageStochasticBlock` also ties the day-ahead declared-dispatch bid
  across the short-period scenarios.

Stochastic scenario sampling uses `pem_extraction`, `scenario_definition` and
`Scen_eps_sampler`. Profiles whose `std` field is missing fall back to a
multiplicative noise `sigma * |mean|`.

## Author

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa
