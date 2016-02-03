#pragma once

#include "Euler.h"
#include <cmath>
#include <armadillo>

struct KelvinHelmholtz : Euler
{
    struct KHSolution : ExactSolution
    {
        double gamma;
        
        double u(double x, double y) const;
        double v(double x, double y) const;
        double rho(double x, double y) const;
        double p(double x, double y) const;
        double rhoE(double x, double y) const;
        
        arma::mat operator()(double x, double y) const;
        
        KHSolution() {};
    };
    
    KelvinHelmholtz(PolyMesh &m, double g);
    
    MeshFn exactSolution(const double t, const int deg);
    MeshFn initialConditions(int deg)
    {
        return exactSolution(0, deg);
    }
};
