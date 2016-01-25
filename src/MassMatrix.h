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
    // Use the product functor to integrate the product of basis functions
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
    
    // Compute the entries of the mass matrix
    MassMatrix(PolyMesh &msh, int deg);
    
    // Compute M x = fn, for given fn, return x
    MeshFn solve(const MeshFn &fn) const;
    // Compute M fn = x for given fn, return x
    MeshFn dot(const MeshFn &fn) const;
};
