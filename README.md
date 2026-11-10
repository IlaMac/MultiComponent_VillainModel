# MultiComponent_VillainModel

This project implements a Monte Carlo simulation for a multi-component Villain model used to study intercomponent Andreev-Bashkin couplings in quantum many-body systems.

The code is written primarily in C++ and is built with CMake. It supports optional MPI, OpenMP, HDF5, Eigen3, spdlog, and test targets.

## Project overview

The repository contains the simulation core, numerical routines, and build configuration needed to compile and run the model:

- `source/` — main C++ implementation and core model logic
- `tests/` — optional test executables configured by CTest
- `cmake/` — project CMake modules and dependency setup
- `CMakeLists.txt` — top-level build configuration
- `build.sh` — convenience script for configuring and building on local machines or cluster environments
- `CMakePresets.json` — preset-based CMake configuration entries
- `conanfile.txt` — Conan dependency declaration

## Requirements

To build the project, you will typically need:

- CMake 3.12 or newer
- A modern C++ compiler (GCC or Clang)
- MPI support (enabled by default in the project)
- OpenMP support
- HDF5, Eigen3, and spdlog libraries
- Optional: h5pp, depending on build configuration

## Quick start

The repository includes a helper script that configures the build in a consistent way:

```bash
./build.sh -h
./build.sh --arch native --build-type Release --make-threads 8 --enable-tests
```

This script supports common options such as compiler selection, debug/release modes, MPI/OpenMP toggles, dependency download method, and target selection.

## Building with CMake directly

You can also configure the project manually:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The project creates an executable named `Villain_2component` as defined in the root CMake configuration.

## Testing

If tests are enabled, CTest is configured through the project build:

```bash
ctest --test-dir build --output-on-failure
```

## Notes

- The build system is designed for both local development and HPC-style clusters.
- The `build.sh` script includes support for machine-specific environment setup and dependency resolution.
- If you are using a cluster environment, the script can be helpful for loading modules before configuring CMake.

## License

This repository does not currently declare a license in the project files visible here. If you plan to distribute or reuse this code, you may want to confirm the repository licensing status with the project owner.
