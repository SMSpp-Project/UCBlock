# UCBlock

A SMS++ module for modeling Unit Commitment problems in electrical power production.

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
[SMS++](https://gitlab.com/smspp/smspp/wikis/custom).

Optionally, install the library in the system with:
```sh
sudo make install
```

### Usage with CMake

After the module is built, you can use it in your CMake project with:
```cmake
find_package(UCBlock)
target_link_libraries(<my_target> SMS++::UCBlock)
```

## Contributing

This section is not ready yet.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

- **Rafael Durbano Lobato**  
  Department of Applied Mathematics  
  State University of Campinas, Brazil

- **Alì Ghezelsoflu**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

- **Niccolò Iardella**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

## License

This section is not ready yet. See SMS++ library for details.

## Disclaimer

The code is currently provided free of charge for academic purposes only.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
