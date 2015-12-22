#pragma once

#include <armadillo>
#include "PolyMesh.h"
#include "BlockMatrix.h"

struct Jacobian : BlockMatrix
{
    const PolyMesh &msh;
    int basisSize;
    int nc;
    
    Jacobian(const PolyMesh &a_msh, int deg, int a_nc);
};
