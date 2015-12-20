#pragma once

#include <armadillo>
#include "Equation.h"

#define kEulerComponents 4

typedef arma::vec::fixed<kEulerComponents> EulerVariables;

struct Euler : Equation
{
    struct FluxDotGradPsi : VecFunctor
    {
        int m;
        int i;
    
        double gamma;
    
        const arma::vec *psi_x, *psi_y;
        const arma::mat *U;
        PolyMesh &msh;
    
        FluxDotGradPsi(PolyMesh &m, double g) : gamma(g), msh(m) {};
    
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
        PolyMesh &msh;
    
        LaxFriedrichsFlux(PolyMesh &m, double g) : gamma(g), msh(m) {};
    
        arma::vec operator()(double x, double y) const;
    };
    
    double gamma;
    int m;
    
    FluxDotGradPsi volumeTerm;
    LaxFriedrichsFlux boundaryTerm;
    
    Euler(PolyMesh &m, double g);
    
    static EulerVariables computeVariables(double x, double y, int m, const arma::mat &U);
    
    static void flux(const EulerVariables &U, double gamma, 
        arma::vec &flux_x, arma::vec &flux_y, double &c, double &u, double &v);
    
    MeshFn assemble(const MeshFn &f);
    
    arma::vec boundaryIntegral(int i, const arma::vec &psi, const MeshFn &U);
    arma::vec volumeIntegral(int i, const arma::vec &psi_x, const arma::vec &psi_y);
};
