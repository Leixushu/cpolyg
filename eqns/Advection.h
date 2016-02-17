#pragma once

#include "Equation.h"

// Two dimensional advection equation, u_t + div(\beta u) = 0
struct Advection : Equation
{
    arma::mat computeVolumeTerm(double x, double y);
    arma::mat computeBoundaryTerm(double x, double y);
    arma::mat computeVolumeJacobian(double x, double y);
    arma::mat computeBoundaryJacobian(double x, double y);
    
    Advection(PolyMesh &a_msh) : Equation(a_msh, 1) {}
};
