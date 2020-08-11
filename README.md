# khover
Computing Khovanov homology and its first derivative.
The project is released under the 2-clause BSD license (see LICENSE).

## Build

### Requirement
- C++17 compatible compiler (maybe C++14 features, generic lambda, and std::optinal will suffice)
- [CMake](http://www.cmake.org)
- [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) matrix library (e.g. libeigen3-dev for Ubuntu).

### Generic instruction
```bash
git clone https://github.com/Junology/khover
cd khover
mkdir build && cd build
cmake ..
make
```
