# replicate_units

Writes a `UCBlock` instance with `K` times the thermal units of a given one:
every unit is copied `K` times, the costs of each copy after the first are
multiplied by a factor drawn uniformly in `[ 1 - p , 1 + p ]`, so that the
copies are not identical, and the demand is multiplied by `K`, so that the
instance keeps the shape of the original one with `K` times the components.

    python replicate_units.py in.nc4 out.nc4 K [ p [ seed ] ]

with `p` defaulting to `0.05` and the seed to `1`; it needs `netCDF4` and
`numpy`. It was written to have instances with many light components, i.e.,
where the subproblems are dynamic programs costing microseconds while the
master problem of a Lagrangian dual grows with their number; on the instances
of the `T-Ramp` family replicated up to `5000` units the partial aggregation
of `BundleSolver` is then worth a factor of `26`.

Only the groups of the thermal units and the demand are touched: every other
group of the file, e.g. `NetworkData`, is copied as it is, so the result is a
single-bus instance if the original one is.
