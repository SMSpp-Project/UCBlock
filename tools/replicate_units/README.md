# replicate_units

Writes a `UCBlock` instance with `K` times the units of a given one: every
unit is copied `K` times, the costs of each copy after the first are
multiplied by a factor drawn uniformly in `[ 1 - p , 1 + p ]`, so that the
copies are not identical, and what counts the units is multiplied by `K` with
them, so that the instance keeps the shape of the original one with `K` times
the components.

    python replicate_units.py in.nc4 out.nc4 K [ p [ seed ] ]

with `p` defaulting to `0.05` and the seed to `1`; it needs `netCDF4` and
`numpy`. It was written to have instances with many light components, i.e.,
where the subproblems are dynamic programs costing microseconds while the
master problem of a Lagrangian dual grows with their number; on the instances
of the `T-Ramp` family replicated up to `5000` units the partial aggregation
of `BundleSolver` is then worth a factor of `26`.

What grows with `K` is of three kinds: the dimensions that count the
replicated components (`NumberUnits`, `NumberElectricalGenerators` and
`NumberStorages`), together with everything indexed by them, which is repeated
`K` times (`GeneratorNode` and the data given per generator, `PollutantRho`
included); the demands, i.e. `ActivePowerDemand`, `ReactivePowerDemand`,
`PrimaryDemand`, `SecondaryDemand` and `InertiaDemand`; and the budgets,
`PollutantBudget` and `PollutantMinBudget`. Leaving any of these as they are
would not give the same instance with `K` times the components, but a
different problem, e.g., one where `K` times the capacity answers the demand
of `n` units, or one where the budget binds `K` times as much.

Everything else is copied as it is, so the result is a single-bus instance if
the original one is, and the units sit on the nodes the original ones sat on.
A file whose other groups, e.g. `NetworkData`, are indexed by one of the
replicated dimensions is refused rather than written half right, as is one
where a quantity would have to be both repeated and scaled.
