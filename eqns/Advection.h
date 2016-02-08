#pragma once

#include "Equation.h"

// Two dimensional advection equation, u_t + div(\beta u) = 0
struct Advection : Equation
{
    struct betaUDotGradPsi : VolumeTermFunctor
    {
        betaUDotGradPsi(PolyMesh &a_msh) : VolumeTermFunctor(a_msh) { n_rows = 1; };
        arma::mat operator()(double x, double y) const;
    };
    
    struct uPsiBetaDotN : NumericalFluxFunctor
    {
        uPsiBetaDotN(PolyMesh &a_msh) : NumericalFluxFunctor(a_msh) { n_rows = 1; };
        arma::mat operator()(double x, double y) const;
    };
    
    struct JacobianBetaUDotGradPsi : VolumeTermJacobianFunctor
    {
        JacobianBetaUDotGradPsi(PolyMesh &a_msh) : VolumeTermJacobianFunctor(a_msh)
        {
            n_rows = 1;
        };
        arma::mat operator()(double x, double y) const;
    };
    
    struct phiPsiBetaDotN : NumericalFluxJacobianFunctor
    {
        phiPsiBetaDotN(PolyMesh &a_msh) : NumericalFluxJacobianFunctor(a_msh) { n_rows = 1; };
        arma::mat operator()(double x, double y) const;
    };
    
    Advection(PolyMesh &m);
    ~Advection();
};
