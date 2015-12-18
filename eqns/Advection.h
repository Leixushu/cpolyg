#pragma once

#include "Equation.h"

struct betaUDotGradPsi : FnFunctor
{
    int m;
    int i;
    
    arma::vec *psi_x, *psi_y;
    arma::vec u;
    PolyMesh &msh;
    
    betaUDotGradPsi(PolyMesh &m) : msh(m) {};
    
    double operator()(double x, double y);
};

struct uPsiBetaDotN : FnFunctor
{
    double nx, ny;
    int m;
    
    int iMinus, iPlus;
    
    arma::vec *psi;
    arma::vec uMinus, uPlus;
    PolyMesh &msh;
    
    uPsiBetaDotN(PolyMesh &m) : msh(m) {};
    
    double operator()(double x, double y);
};

struct Advection : Equation
{
    double beta_x, beta_y;
    betaUDotGradPsi volumeTerm;
    uPsiBetaDotN boundaryTerm;
    
    Advection(PolyMesh & m);
    
    MeshFn assemble(const MeshFn &f);
    
    double volumeIntegral(int i, arma::vec &psi_x, arma::vec &psi_y);
    double boundaryIntegral(int i, arma::vec &psi, const MeshFn &u);
};
