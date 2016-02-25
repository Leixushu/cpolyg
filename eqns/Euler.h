#pragma once

#include <armadillo>
#include "Equation.h"

#define kEulerComponents 4

typedef arma::vec::fixed<kEulerComponents> EulerVariables;
typedef arma::mat::fixed<kEulerComponents, kEulerComponents> EulerJacobian;

struct Euler : Equation
{
    arma::mat fluxFunction(const arma::vec &vars, double x, double y);
    arma::vec numericalFluxFunction(const arma::vec &varsMinus, 
        const arma::vec &varsPlus, 
        double x, double y, double nx, double ny);
    arma::cube fluxJacobian(const arma::vec &vars, double x, double y);
    arma::mat numericalFluxJacobian(const arma::vec &varsMinus, 
        const arma::vec &varsPlus, 
        double x, double y, double nx, double ny, int sgn);
    
    Euler(PolyMesh &m, BoundaryConditions a_bc, double g)
    : Equation(m, a_bc, kEulerComponents), gamma(g) { }
    
    EulerJacobian Id = arma::eye(kEulerComponents, kEulerComponents);
    double gamma;
    
protected:
    
    void flux(const EulerVariables &U, arma::vec &flux_x, arma::vec &flux_y, 
              double &c, double &u, double &v);
    EulerVariables alphaDerivative(const EulerVariables &vars1, 
                                   const EulerVariables &vars2,
                                   double &alpha);
};
