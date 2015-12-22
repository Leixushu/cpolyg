#pragma once

#include <armadillo>
#include <string>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Legendre.h"
#include "BlockMatrix.h"

struct MassMatrix : BlockMatrix
{
    struct ProductFunctor : FnFunctor
    {
        arma::vec phi, psi;
        int i;
        int m;
        PolyMesh &msh;
    
        ProductFunctor(PolyMesh &a_msh) : msh(a_msh) {};
    
        double operator()(double x, double y) const
        {
            double xx, yy;
            msh.getLocalCoordinates(i, x, y, xx, yy);
        
            return Leg2D(xx, yy, m, phi)*Leg2D(xx, yy, m, psi);
        };
    };
    
    PolyMesh  &msh;
    int deg;
    int basisSize;
    
    MassMatrix(PolyMesh &msh, int deg);
    
    MeshFn solve(const MeshFn &fn) const;
};
