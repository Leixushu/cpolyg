#pragma once

#include <armadillo>
#include "Equation.h"

#define kEulerComponents 4

typedef arma::vec::fixed<kEulerComponents> EulerVariables;
typedef arma::mat::fixed<kEulerComponents, kEulerComponents> EulerJacobian;

struct Euler : Equation
{
    arma::mat computeVolumeTerm(double x, double y);
    arma::mat computeBoundaryTerm(double x, double y);
    arma::mat computeVolumeJacobian(double x, double y);
    arma::mat computeBoundaryJacobian(double x, double y);
    
    Euler(PolyMesh &m, BoundaryConditions a_bc, double g);
    
    
    EulerJacobian Id = arma::eye(kEulerComponents, kEulerComponents);
    double gamma;
    
    
protected:

    static EulerVariables computeVariables(double x, double y, 
                                           int m, const arma::mat &U);
    
    void flux(const EulerVariables &U, arma::vec &flux_x, arma::vec &flux_y, 
              double &c, double &u, double &v);
    void fluxJacobian(const EulerVariables &U, EulerJacobian &J1, EulerJacobian &J2);
    EulerVariables alphaDerivative(const EulerVariables &vars1, 
                                   const EulerVariables &vars2,
                                   double &alpha);
};
