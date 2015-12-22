#pragma once

#include <cstdlib>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "BlockMatrix.h"

struct Equation
{
    PolyMesh &msh;
    
    Equation(PolyMesh &m) : msh(m) { };
    
    virtual MeshFn assemble(const MeshFn &f, double t) = 0;
    virtual BlockMatrix jacobian(const MeshFn &f, double t)
    {
        abort();
        return BlockMatrix();
    };
    
    virtual ~Equation() {};
};
