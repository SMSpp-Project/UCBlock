# Data Converter Tool

This is a tool to convert UC data in the Energy Communities setting from
csv to NetCDF.

## Getting started

The main script file is `csv2nc4.jl` that, if run without specific options, 
i.e.:

```
julia csv2nc4.jl
```

will produce a nc4 file optimized wrt the redundancy of the input data.
In another sense, if all the `NetworkBlock`s of the problem conceptually 
contain the same data about the market (i.e., sell, buy and reward price), 
the script does not will physically generate it all, but it will store
the data just one time in the father block, i.e., `Block_0`.
In this case, the `ActiveDemand` will be stored in the father block, i.e., 
`Block_0`, or in its respectively `NetworkBlock` otherwise.

If you want to force the physical creation of all the `NetworkBlock`s needed 
to represent the problem, e.g., for testing reasons, you can use the option
`-with-network-blocks`, i.e.:

```
julia csv2nc4.jl -with-network-blocks
```