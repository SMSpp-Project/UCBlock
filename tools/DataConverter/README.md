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
contain the same data about the market (i.e., sell price and buy price), 
the script does not will physically generate it all, but it will store
the data just one time in the father block, i.e., `Block_0`.
In this case, the `ActiveDemand` and the `ConstTerm` will be stored in the 
father block, i.e., `Block_0`, as `ActivePowerDemand` and `ConstTerm` 
respectively.

If you want to force the physical creation of all the `NetworkBlock`s needed 
to represent the problem, e.g., for testing reasons, you can use the option
`-with-network-blocks`, i.e.:

```
julia csv2nc4.jl -with-network-blocks
```

In this case, by default, the `ActiveDemand` and the `ConstTerm` will be 
stored in the `NetworkBlock` to which they refer.

TODO : consider to add other two options, i.e., `-store-demand-in-father`
and `-store-const-term-in-father`, that case of `-with-network-blocks`, 
enforce the storage of the `ActiveDemand` and the `ConstTerm` respectively just 
one time in the father block, i.e., `Block_0`, always as `ActivePowerDemand` 
and `ConstTerm` respectively.