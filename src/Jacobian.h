#pragma once

#include <armadillo>
#include "PolyMesh.h"
#include "BlockMatrix.h"
#include "MeshFn.h"
#include "MassMatrix.h"

struct Jacobian : BlockMatrix
{
    PolyMesh &msh;
    int nc;
    int deg;
    
    Jacobian(PolyMesh &a_msh) : msh(a_msh) {};
    Jacobian(PolyMesh &a_msh, int a_deg, int a_nc);
    
    MeshFn dot(const MeshFn &x);
    MeshFn solve(const MeshFn &b, Preconditioner &pc);
    
    Jacobian& operator +=(MassMatrix &M);
};
