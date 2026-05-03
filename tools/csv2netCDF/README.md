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
- `energy_community_model_NA`     (i.e., No Asset case)
- `energy_community_model_NC`     (i.e., No Cooperative case)

The same names with the `_sto` suffix enable scenario sampling (the
deterministic UCBlock is still produced; an additional `TSSB_*.nc4` is
written that references it via the `filename` attribute).

It is also possible to choose whether to include the generator thermal using
the appropriate flag, i.e.:

```sh
julia csv2nc4.jl [yml] --with-thermal-blocks
```

Finally, for tests purposes, it can take an additional parameter to enforce
the generation of the physical ECNetworkBlock(s), i.e.;

```sh
julia csv2nc4.jl [yml] --with-network-blocks
```

The convenience wrapper `gen-all-nc4` runs the deterministic and stochastic
batches over CO, NA and NC, with all four flag combinations.

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
