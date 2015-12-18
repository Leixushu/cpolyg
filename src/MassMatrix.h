#pragma once

#include <armadillo>
#include <string>
#include "PolyMesh.h"
#include "MeshFn.h"

struct MassMatrix
{
    arma::sp_mat matrix;
    PolyMesh  &msh;
    int deg;
    
    MassMatrix(PolyMesh &msh, int deg);
    
    MeshFn solve(const MeshFn &fn) const;
    void spy(std::string filename);
};
