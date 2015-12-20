#pragma once

#include <armadillo>
#include <string>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Legendre.h"

struct MassMatrix
{
    struct ProductFunctor : FnFunctor
    {
        arma::vec phi, psi;
        int i;
        int m;
        PolyMesh &msh;
    
        ProductFunctor(PolyMesh &m) : msh(m) {};
    
        double operator()(double x, double y) const
        {
            double xx, yy;
            msh.getLocalCoordinates(i, x, y, xx, yy);
        
            return Leg2D(xx, yy, m, phi)*Leg2D(xx, yy, m, psi);
        };
    };
    
    arma::sp_mat matrix;
    PolyMesh  &msh;
    int deg;
    
    MassMatrix(PolyMesh &msh, int deg);
    
    MeshFn solve(const MeshFn &fn) const;
    void spy(std::string filename);
};
