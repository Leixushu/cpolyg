#pragma once

#include "Euler.h"
#include <cmath>
#include <armadillo>

struct EulerVortex : Euler
{
    struct VortexSolution : ExactSolution
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
        
        arma::vec operator()(double x, double y) const;
        
        VortexSolution() {};
    };
    
    EulerVortex(PolyMesh &m, double g);
    
    MeshFn exactSolution(const double t, const int deg);
    MeshFn initialConditions(int deg)
    {
        return exactSolution(0, deg);
    }
};
