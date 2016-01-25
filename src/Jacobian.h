#pragma once

#include <armadillo>
#include "PolyMesh.h"
#include "BlockMatrix.h"
#include "MeshFn.h"
#include "MassMatrix.h"

enum Solver
{
    kGMRESSolver,
    kJacobiSolver
};

struct Jacobian : BlockMatrix
{
    PolyMesh &msh;
    int nc;
    int deg;
    
    Jacobian(PolyMesh &a_msh) : msh(a_msh) {};
    Jacobian(PolyMesh &a_msh, int a_deg, int a_nc);
    
    MeshFn dot(const MeshFn &x);
    MeshFn solve(const MeshFn &b, Preconditioner &pc, Solver s=kGMRESSolver);
    
    Jacobian& operator +=(const MassMatrix &M);
    Jacobian& operator=(const Jacobian &J2);
};
