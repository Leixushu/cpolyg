#pragma once

#include <armadillo>
#include <string>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Legendre.h"
#include "BlockMatrix.h"

// Implement the mass matrix as a special type of block matrix
struct MassMatrix : BlockMatrix
{
    std::vector<arma::mat> LU;
    std::vector<arma::Col<int>> ipvt;
    
    PolyMesh  &msh;
    int deg;
    int nc;
    
    // Compute the entries of the mass matrix
    MassMatrix(PolyMesh &a_msh, int a_deg, int a_nc);
    
    // Compute M x = fn, for given fn, return x
    MeshFn solve(const MeshFn &fn);
    // Compute M fn = x for given fn, return x
    MeshFn matvec(const MeshFn &fn) const;
};
