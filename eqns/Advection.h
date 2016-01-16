#pragma once

#include "Equation.h"

/** \addtogroup Equations
    @{
*/

/// Two dimensional advection equation, \f$ u_t + \nabla\cdot(\beta u) = 0 \f$
struct Advection : Equation
{
    struct betaUDotGradPsi : VolumeTermFunctor
    {
        betaUDotGradPsi(PolyMesh &a_msh) : VolumeTermFunctor(a_msh) { nc = 1; };
        arma::mat operator()(double x, double y) const;
    };
    
    struct uPsiBetaDotN : NumericalFluxFunctor
    {
        uPsiBetaDotN(PolyMesh &a_msh) : NumericalFluxFunctor(a_msh) { nc = 1; };
        arma::mat operator()(double x, double y) const;
    };
    
    struct JacobianBetaUDotGradPsi : VolumeTermJacobianFunctor
    {
        JacobianBetaUDotGradPsi(PolyMesh &a_msh) : VolumeTermJacobianFunctor(a_msh)
        {
            nc = 1;
        };
        arma::mat operator()(double x, double y) const;
    };
    
    struct phiPsiBetaDotN : NumericalFluxJacobianFunctor
    {
        phiPsiBetaDotN(PolyMesh &a_msh) : NumericalFluxJacobianFunctor(a_msh) { nc = 1; };
        arma::mat operator()(double x, double y) const;
    };
    
    Advection(PolyMesh &m);
    ~Advection();
};

/**@}*/
