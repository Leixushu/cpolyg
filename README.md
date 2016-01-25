# README #

C++ code for polygonal discontinuous Galerkin method.  
Will Pazner  
December 20, 2015

The code uses:

- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html), from Jonathan Shewchuk
- [Voro++](http://math.lbl.gov/voro++/), from Chris Rycroft
- GMRES code from [Per-Olof Persson](http://persson.berkeley.edu)

## Build Instructions

`cpolyg` can be built by invoking `make` in the main directory. Edit the Makefile to 
specify the compiler, (`CC` and `CXX`), and the system include and library directories,
which are required to link against Armadillo, GSL, SuperLU, BLAS, and LAPACK.
