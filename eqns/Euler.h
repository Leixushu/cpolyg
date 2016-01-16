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
        
        virtual double f(double x, double y) const = 0;
        virtual double u(double x, double y) const = 0;
        virtual double v(double x, double y) const = 0;
        virtual double rho(double x, double y) const = 0;
        virtual double p(double x, double y) const = 0;
        virtual double rhoE(double x, double y) const = 0;
        
        virtual arma::mat operator()(double x, double y) const = 0;
        
        ExactSolution() { nc = kEulerComponents; t = 0; };
        virtual ~ExactSolution() { };
    };
    
    struct FluxDotGradPsi : VolumeTermFunctor
    {
        double gamma;
    
        FluxDotGradPsi(PolyMesh &a_msh, double a_gamma) : VolumeTermFunctor(a_msh)
        {
            gamma = a_gamma;
            nc = kEulerComponents;
        }
    
        arma::mat operator()(double x, double y) const;
    };

    struct LaxFriedrichsFlux : NumericalFluxFunctor
    {
        double gamma;
        ExactSolution *exact;
    
        LaxFriedrichsFlux(PolyMesh &a_msh, double a_gamma) : NumericalFluxFunctor(a_msh)
        {
            gamma = a_gamma;
            nc = kEulerComponents;
            exact = NULL;
        }
    
        arma::mat operator()(double x, double y) const;
    };
    
    ExactSolution *exact;
    
    double gamma;
    
    static EulerVariables computeVariables(double x, double y, int m, const arma::mat &U);
    static void flux(const EulerVariables &U, double gamma, 
        arma::vec &flux_x, arma::vec &flux_y, double &c, double &u, double &v);
    
    Euler(PolyMesh &m, double g);
    ~Euler();
};
