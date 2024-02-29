# UCBlock

A SMS++ module for modeling Unit Commitment (UC) problems in electrical power
production.

The core of the module is the `UCBlock` class which represents the UC problem.
The formulation is very flexible in that `UCBlock` only knows that is has a
bunch of generating units, each one a concrete class deriving from the abstract
base class `UnitBlock`; several of these are available, such as
`ThermalUnitBlock`, `HydroUnitBlock`, `BatteryUnitBlock` and others. Also,
`UCBlock` knows that energy must flow between generating units and consumption
points through an energy network, represented by concrete class deriving from
the abstract base class `NetworkBlock` (unless there is no networ, i.e., the
"bus" case, which is handled directly by `UCBlock`); some of these are
available, such as `DCNetworkBlock` for the linear DC or HVDC (or hybrid)
cases and `ECNetworkBlock` for Energy Communities having to share the
energy between users and then with the external grid. Other kinds of
units and networks can easily be added, and specialised solution methods for
certaint units and networks (e.g., `ThermalUnitDPSolver` for
`ThermalUnitBlock`) can be developed.


## Getting started

These instructions will let you build UCBlock on your system.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
sudo make install
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(UCBlock)
target_link_libraries(<my_target> SMS++::UCBlock)
```

### Running the tests with CMake

Some unit tests will be built with the library. Launch `ctest` from the
build directory to run them. To disable them, set the option
`BUILD_TESTING` to `OFF`.

> **Note:**
> CMake will fetch and build it automatically.

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `UCBlock` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the [tester for
ThermalUnitDPSolver](https://gitlab.com/smspp/tests/-/tree/develop/ThermalUnitBlock_Solver?ref_type=heads). The makefiles in turn recursively include all the
required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at their
default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Tools

We provide some tool to generate input data for UCBlock:

- [a converter from text-based formats to netCDF](tools/nc4generator.cpp)
  that can be used to produce netCDF versions of the instances produced by
  [classical random generators](https://commalab.di.unipi.it/datasets/UC)

- [a Matlab-based data generator](tools/DataGenerator/README.md)

- [a converter from .yml and .csv data files](tools/DataConverter/README.md)
  used to describe UC instances corresponding to Energy Community design
  problems used in the [EnergyCommunity.jl JuMP
  package](https://github.com/SPSUnipi/EnergyCommunity.jl)



## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/ucblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

- **Tiziano Bacci**  
  Istituto di Analisi dei Sistemi ed Informatica "A. Ruberti"  
  Consiglio Nazionale delle Ricerche

### Previous Contributors

- **Ali Ghezelsoflu**  
  Dipartimento di Informatica  
  Università di Pisa

- **Niccolo' Iardella**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
