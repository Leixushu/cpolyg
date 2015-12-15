#pragma once

#include "Equation.h"

struct Advection : Equation
{
    MeshFn assemble(MeshFn &f);
};
