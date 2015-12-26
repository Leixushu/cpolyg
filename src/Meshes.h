#pragma once

#include "PolyMesh.h"

PolyMesh quadRectangle(double h, double width, double height);
PolyMesh hexRectangle(double h, double width, double height);
PolyMesh triRectangle(double h, double width, double height);
PolyMesh honeycombRectangle(double h, double width, double height);

PolyMesh honeycombUnitSquare(double h);
PolyMesh quadUnitSquare(double h);
PolyMesh hexUnitSquare(double h);
PolyMesh triUnitSquare(double h);
