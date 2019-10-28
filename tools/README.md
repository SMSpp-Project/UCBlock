# UCBlock Tools

A set of tools and examples that use UCBlock module.
At the moment we provide only a single Thermal Unit solver.

## Getting started

These instructions will let you build UCBlock Tools on your system.

### Requirements

- UCBlock

At the moment, the Thermal Unit solver also requires:

- MILPSolver

### Build and install

Configure and build with:
```sh
mkdir build
cd build
cmake ..
make
```

## Usage

### Thermal Unit solver

```sh
Usage: thermalunit_solver [options] <nc4-file>

Options:
  -s <solver>, --solver <solver>  Choose solver.
                                  Available solvers are: cplex, dp.
  -w <file>, --writelp <file>     Write LP problem on file.
  -n <file>, --nc4problem <file>  Write nc4 problem on file.
  -h, --help                      Print this help.
```

The input netCDF file must be a block file.
At the moment, only `cplex` solver is available.

## Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

- **Niccolò Iardella**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

## License

See UCBlock for details.

## Disclaimer

The code is currently provided free of charge for academic purposes only.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
