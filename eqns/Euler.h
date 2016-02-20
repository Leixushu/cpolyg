#pragma once

#include <armadillo>
#include "Equation.h"

#define kEulerComponents 4

typedef arma::vec::fixed<kEulerComponents> EulerVariables;
typedef arma::mat::fixed<kEulerComponents, kEulerComponents> EulerJacobian;

struct Euler : Equation
{
    struct ExactSolution : VecFunctor
    {
        double t;
        
        virtual double u(double x, double y) const = 0;
        virtual double v(double x, double y) const = 0;
        virtual double rho(double x, double y) const = 0;
        virtual double p(double x, double y) const = 0;
        virtual double rhoE(double x, double y) const = 0;
        
        virtual arma::mat operator()(double x, double y) const = 0;
        
        ExactSolution() : VecFunctor(kEulerComponents) { t = 0; };
        virtual ~ExactSolution() { };
    };
    
    ExactSolution *exact = nullptr;
    EulerJacobian Id = arma::eye(kEulerComponents, kEulerComponents);
    double gamma;
    
    Euler(PolyMesh &m, double g) : Equation(m, kEulerComponents), gamma(g) { };
    ~Euler();
    
    arma::mat computeVolumeTerm(double x, double y);
    arma::mat computeBoundaryTerm(double x, double y);
    arma::mat computeVolumeJacobian(double x, double y);
    arma::mat computeBoundaryJacobian(double x, double y);
    
    static EulerVariables computeVariables(double x, double y, int m, const arma::mat &U);
    
    void flux(const EulerVariables &U, arma::vec &flux_x, arma::vec &flux_y, 
              double &c, double &u, double &v);
    void fluxJacobian(const EulerVariables &U, EulerJacobian &J1, EulerJacobian &J2);
    EulerVariables alphaDerivative(const EulerVariables &vars1, 
                                   const EulerVariables &vars2,
                                   double &alpha);
    
    Jacobian jacobian(const MeshFn &f, double t);
    MeshFn assemble(const MeshFn &f, double t);
};
