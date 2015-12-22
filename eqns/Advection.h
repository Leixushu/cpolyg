#pragma once

#include "Equation.h"

struct Advection : Equation
{
    struct betaUDotGradPsi : FnFunctor
    {
        int m;
        int i;
    
        arma::vec *psi_x, *psi_y;
        arma::vec u;
        PolyMesh &msh;
    
        betaUDotGradPsi(PolyMesh &a_msh) : msh(a_msh) {};
        
        double operator()(double x, double y) const;
    };

    struct uPsiBetaDotN : FnFunctor
    {
        double nx, ny;
        int m;
    
        int iMinus, iPlus;
    
        arma::vec *psi;
        arma::vec uMinus, uPlus;
        PolyMesh &msh;
    
        uPsiBetaDotN(PolyMesh &a_msh) : msh(a_msh) {};
    
        double operator()(double x, double y) const;
    };
    
    betaUDotGradPsi volumeTerm;
    uPsiBetaDotN boundaryTerm;
    
    Advection(PolyMesh &m);
    
    MeshFn assemble(const MeshFn &f, double t = 0);
    
    
    
    double volumeIntegral(int i, arma::vec &psi_x, arma::vec &psi_y);
    double boundaryIntegral(int i, arma::vec &psi, const MeshFn &u);
};
