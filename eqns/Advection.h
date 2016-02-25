#pragma once

#include "Equation.h"

// Two dimensional advection equation, u_t + div(\beta u) = 0
struct Advection : Equation
{
    arma::mat fluxFunction(const arma::vec &vars, double x, double y);
    arma::vec numericalFluxFunction(const arma::vec &varsMinus, 
        const arma::vec &varsPlus, 
        double x, double y, double nx, double ny);
    arma::cube fluxJacobian(const arma::vec &vars, double x, double y);
    arma::mat numericalFluxJacobian(const arma::vec &varsMinus,
        const arma::vec &varsPlus, 
        double x, double y, double nx, double ny, int sgn);
    
    virtual double betaX(double x, double y) const;
    virtual double betaY(double x, double y) const;
    
    Advection(PolyMesh &a_msh, BoundaryConditions a_bc)
    : Equation(a_msh, a_bc, 1) {}
};

struct ConstantAdvection : Advection
{
    double bx, by;
    
    ConstantAdvection(PolyMesh &a_msh, BoundaryConditions a_bc, 
                      double a_bx, double a_by)
    : Advection(a_msh, a_bc), bx(a_bx), by(a_by)
    { }
    
    double betaX(double x, double y) const
    {
        return bx;
    }
    
    double betaY(double x, double y) const
    {
        return by;
    }
};
