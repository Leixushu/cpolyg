#pragma once

#include "PolyMesh.h"

// Routines for generating commonly-used meshes

// tile (width)x(height) rectangle by quadrilaterals, hexagons, or triangles
PolyMesh quadRectangle(double h, double width, double height);
PolyMesh hexRectangle(double h, double width, double height);
PolyMesh triRectangle(double h, double width, double height);
PolyMesh honeycombRectangle(double h, double width, double height);

// convenience routines for the unit square
PolyMesh honeycombUnitSquare(double h);
PolyMesh quadUnitSquare(double h);
PolyMesh hexUnitSquare(double h);
PolyMesh triUnitSquare(double h);
