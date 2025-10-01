# 2D Nonlinear Poisson Solver

A C++ implementation for solving 2D nonlinear Poisson equations using the finite element method with Newton iteration.

## Features

- **Finite Element Method**: Support for triangular and quadrilateral elements
- **Nonlinear Solver**: Newton iteration with adaptive tolerance
- **Sparse Matrices**: Efficient memory usage for large problems
- **Multiple Element Types**: Linear triangles and quads
- **VTK Output**: Visualization with ParaView
- **Comprehensive Testing**: Unit tests and analytical verification
- **Modern C++**: CMake build system, Google Test, Doxygen docs

## Building

```bash
# Clone and build
git clone <repository>
cd Poisson2DSolver
mkdir build && cd build
cmake .. -DUSE_SPARSE_MATRIX=ON -DBUILD_TESTS=ON
make -j4

# Run tests
./Poisson2D_tests

# Run solver
./Poisson2DSolver ../examples/input.json