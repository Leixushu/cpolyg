#pragma once

#include "PolyMesh.h"

// Routines for generating commonly-used meshes

// tile (width)x(height) rectangle by quadrilaterals, hexagons, or triangles
PolyMesh quadRectangle(double h, double width, double height);
PolyMesh perturbedQuadRectangle(double h, double p, double width, double height);
PolyMesh perturbedTriRectangle(double h, double p, double width, double height);
PolyMesh hexRectangle(double h, double width, double height);
PolyMesh triRectangle(double h, double width, double height);
PolyMesh honeycombRectangle(double h, double width, double height);
PolyMesh periodicRectangle(int Nx, int Ny, double width, double height);

PolyMesh rectangle1D(double h, double width, double height);

// convenience routines for the unit square
PolyMesh honeycombUnitSquare(double h);
PolyMesh quadUnitSquare(double h);
PolyMesh hexUnitSquare(double h);
PolyMesh triUnitSquare(double h);
