#pragma once

#include "Equation.h"

struct betaUDotGradPsi : FnFunctor
{
    double beta_x, beta_y;
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
    double beta_x, beta_y, nx, ny;
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
    
    Advection(PolyMesh & m)
    : Equation(m), volumeTerm(betaUDotGradPsi(m)), boundaryTerm(uPsiBetaDotN(m))
    {
        beta_x = 1;
        beta_y = 0;
    };
    
    MeshFn assemble(MeshFn &f);
    
    double volumeIntegral(int i, arma::vec &psi_x, arma::vec &psi_y);
    double boundaryIntegral(int i, arma::vec &psi, MeshFn &u);
};
