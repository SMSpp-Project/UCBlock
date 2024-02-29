# Data Converter Tool

This is a tool to convert UC data in the Energy Communities setting from
csv to NetCDF.

## Getting started

The main script file is `csv2nc4.jl` that optionally takes in input the 
following parameters:

```sh
julia csv2nc4.jl [yml]
```

where `yml` can be one of the followings:

- `energy_community_model_CO` (i.e., Cooperative case; the default, if none is given)
- `energy_community_model_NA` (i.e., No Asset case)
- `energy_community_model_NC` (i.e., No Cooperative case)

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

## Author

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa