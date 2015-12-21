#pragma once

#include <armadillo>
#include "Equation.h"

#define kEulerComponents 4

typedef arma::vec::fixed<kEulerComponents> EulerVariables;

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
        
        virtual arma::vec operator()(double x, double y) const = 0;
        
        ExactSolution() { nc = kEulerComponents; };
        virtual ~ExactSolution() { };
    };
    
    struct FluxDotGradPsi : VecFunctor
    {
        int m;
        int i;
    
        double gamma;
    
        const arma::vec *psi_x, *psi_y;
        const arma::mat *U;
        PolyMesh &msh;
    
        FluxDotGradPsi(PolyMesh &a_msh, double a_gamma) : gamma(a_gamma), msh(a_msh)
        {
            nc = kEulerComponents;
        }
    
        arma::vec operator()(double x, double y) const;
    };

    struct LaxFriedrichsFlux : VecFunctor
    {
        double nx, ny;
        int m;
        
        double gamma;
        
        int iMinus, iPlus;
    
        const arma::vec *psi;
        const arma::mat *UMinus, *UPlus;
        
        ExactSolution *exact;
        
        PolyMesh &msh;
    
        LaxFriedrichsFlux(PolyMesh &a_msh, double a_gamma) : gamma(a_gamma), msh(a_msh)
        {
            nc = kEulerComponents;
        }
    
        arma::vec operator()(double x, double y) const;
    };
    
    ExactSolution *exact;
    
    double gamma;
    int m;
    
    FluxDotGradPsi volumeTerm;
    LaxFriedrichsFlux boundaryTerm;
    
    Euler(PolyMesh &m, double g);
    ~Euler() { if(exact) delete exact; }
    
    static EulerVariables computeVariables(double x, double y, int m, const arma::mat &U);
    
    static void flux(const EulerVariables &U, double gamma, 
        arma::vec &flux_x, arma::vec &flux_y, double &c, double &u, double &v);
    
    MeshFn assemble(const MeshFn &f, double t);
    
    arma::vec boundaryIntegral(int i, const arma::vec &psi, const MeshFn &U);
    arma::vec volumeIntegral(int i, const arma::vec &psi_x, const arma::vec &psi_y);
};
