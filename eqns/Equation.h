#pragma once

#include "PolyMesh.h"
#include "MeshFn.h"

struct Equation
{
    PolyMesh &msh;
    
    Equation(PolyMesh &m) : msh(m) { };
    
    virtual MeshFn assemble(const MeshFn &f) = 0;
};
