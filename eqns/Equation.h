#pragma once

#include <cstdlib>
#include <armadillo>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Jacobian.h"

struct Equation
{
    typedef arma::mat (Equation::*Integrand)(double x, double y);
    
    struct IntegrandFunctor : VecFunctor
    {
        Equation &eqn;
        Integrand integ;
        
        IntegrandFunctor(Equation &a_eqn, Integrand a_integ, int a_n_rows, int a_n_cols)
        : VecFunctor(a_n_rows, a_n_cols), eqn(a_eqn), integ(a_integ) { }
        
        arma::mat operator()(double x, double y) const
        {
            return (eqn.*integ)(x, y);
        }
    };
    
    PolyMesh &msh;
    
    IntegrandFunctor *volumeTerm;
    IntegrandFunctor *boundaryTerm;
    
    IntegrandFunctor *volumeJacobian;
    IntegrandFunctor *boundaryJacobian;
    
    int nc;
    
    // needed during assemble/jacobian routines
    arma::vec phi, psi, psi_x, psi_y;
    arma::mat UMinus, UPlus, U, UNeighbor;
    int m, iPlus, iMinus, iPhi, iPsi, neighbor;
    double nx, ny;
    
    
    Equation(PolyMesh &a_msh, int a_nc);
    virtual ~Equation();
    
    arma::vec boundaryIntegral(int i, const MeshFn &u, int deg);
    arma::vec volumeIntegral(int i, int deg);
    
    virtual MeshFn assemble(const MeshFn &f, double a_t);
    virtual Jacobian jacobian(const MeshFn &f, double a_t);
    
    virtual arma::mat computeVolumeTerm(double x, double y) = 0;
    virtual arma::mat computeBoundaryTerm(double x, double y) = 0;
    virtual arma::mat computeVolumeJacobian(double x, double y) = 0;
    virtual arma::mat computeBoundaryJacobian(double x, double y) = 0;
};
