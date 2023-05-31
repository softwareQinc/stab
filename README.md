# stab
Efficient simulation of Clifford circuits using the method described in 
[arXiv:2109.08629 [quant-ph]](https://arxiv.org/abs/2109.08629)

## Pre-requisites

- C++17 compliant compiler, e.g., [GNU gcc](https://gcc.gnu.org/)
  (`sudo apt install build-essential` to install on Ubuntu/Debian Linux)
- [CMake](https://cmake.org/) build system
  (`sudo apt install cmake` to install on Ubuntu/Debian Linux, or `brew install cmake` to install on macOS)
- [Eigen3](https://eigen.tuxfamily.org/index.php) matrix library version 3.4 or later
  (`sudo apt install libeigen3-dev` to install on Ubuntu/Debian Linux, or `brew install eigen` to install on macOS)
- [Quantum++](https://github.com/softwareqinc/qpp) quantum computing library
  (`git clone https://github.com/softwareqinc/qpp && cd qpp && cmake -B build && sudo cmake --build build --target install`
  to install on all platforms, or `brew install quantum++` to install on macOS)

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
