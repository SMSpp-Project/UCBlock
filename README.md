# UCBlock
> TODO

## Requirements

The project requires the [core SMS++ library](https://gitlab.com/frangio68/sms_plus_plus).
It has the same requirements of the core SMS++ project (Boost, Eigen3 and netCDF).

## Build and install

- Clone the project from the repository and navigate inside its main directory.

- If you installed the requirements you should be fine. Configure the project with:
```sh
mkdir build
cd build
cmake ..
```

  If you can't or wont install the required libraries, you will need to specify
  their custom path in the [`CMakeCustom.txt`](CMakeCustom.txt), e.g.:
```cmake
...
# Boost main directory
set(BOOST_ROOT  "/my/custom/path/to/boost/")
...
```

  See [`CMakeCustom.txt`](CMakeCustom.txt) for other ways to customize the configuration of SMS++ projects.

- You can now build the library:
```sh
make
```

- Optionally, you can install the library with:
```sh
sudo make install
```

- After the library is configured and built, you can use it in your CMake project with:
```cmake
find_package(MILPSolver)
target_link_libraries(<my_target> SMS++::MILPSolver)
```

### Run Google Tests

If `enable_testing()` is not commented in the customization file, some simple Google Tests will
be built with the library. You can run them with the following command from the `build` directory:
```sh
ctest
```

## Tools

- `thermalunit_solver`: a tool for solving simple 1-UC problems passed as netCDF files.
```
Usage: thermalunit_solver [options] <file>

-s <solver>, --solver <solver>  Choose solver.
                                Available solvers are: cplex, dp.
-w <file>, --writelp <file>     Write LP problem on file.
-h, --help                      Print this help.
```
  
