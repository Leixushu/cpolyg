#pragma once

#include "Euler.h"
#include <cmath>
#include <armadillo>

struct VortexSolution : VecFunctor
{
    double gamma;
    double x0, y0;
    double theta, epsilon, r_c, MInf, uInf, rhoInf, uBar, vBar, pInf;

    double f(double x, double y) const;
    double u(double x, double y) const;
    double v(double x, double y) const;
    double rho(double x, double y) const;
    double p(double x, double y) const;
    double rhoE(double x, double y) const;
    
    arma::mat operator()(double x, double y) const;
    
    VortexSolution(double a_gamma);
    
    MeshFn interpolated(const PolyMesh &msh, const int deg);
};
