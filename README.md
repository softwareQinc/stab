# stab
Efficient simulation of Clifford circuits using the method described in Niel de Beaudrap and Steven Herbert's paper ["Fast Stabiliser Simulation with Quadratic Form Expansions"](https://quantum-journal.org/papers/q-2022-09-15-803/). Most functions are taken directly from the pseudocode described in that paper.

## Pre-requisites

- C++17 compliant compiler, e.g., [GNU gcc](https://gcc.gnu.org/)
  (`sudo apt install build-essential` to install on Ubuntu/Debian Linux)
- [CMake](https://cmake.org/) build system
  (`sudo apt install cmake` to install on Ubuntu/Debian Linux, or `brew install cmake` to install on macOS)
- [Eigen3](https://eigen.tuxfamily.org/index.php) matrix library version 3.4 or later
  (`sudo apt install libeigen3-dev` to install on Ubuntu/Debian Linux, or `brew install eigen` to install on macOS)

## Optional

- [Quantum++](https://github.com/softwareqinc/qpp) quantum computing library

Execute

```bash
git clone https://github.com/softwareqinc/qpp 
cd qpp 
cmake -B build -DUSE_QPP=ON 
sudo cmake --build build --target install
```

to install on all platforms, or `brew install quantum++` to install on macOS, or `sudo pkg install quantum++` to 
install on FreeBSD. 

## Setup

From inside the project's root directory, execute

```bash
cmake -B build
cmake --build build --parallel 8
```

## Unit testing

To build the unit tests, execute

```bash
cmake --build build --target unit_tests --parallel 8
```

To run the unit tests, execute

```bash
ctest --test-dir build
```
