#pragma once
//#include <cstdlib>
#include <armadillo>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Jacobian.h"
#include "BoundaryConditions.h"

struct Equation
{
    PolyMesh &msh;
    
    struct IntegrandFunctor;
    IntegrandFunctor *volumeTerm;
    IntegrandFunctor *boundaryTerm;
    IntegrandFunctor *volumeJacobian;
    IntegrandFunctor *boundaryJacobian;
    
    BoundaryConditions bc;
    
    // number of components (e.g. 1 for advection, 4 for 2D Euler)
    int nc;
    
    Equation(PolyMesh &a_msh, BoundaryConditions a_bc, int a_nc);
    virtual ~Equation();
    
    virtual MeshFn assemble(const MeshFn &f, double a_t);
    virtual Jacobian jacobian(const MeshFn &f, double a_t);
    
protected:
    
    virtual arma::mat computeVolumeTerm(double x, double y) = 0;
    virtual arma::mat computeBoundaryTerm(double x, double y) = 0;
    virtual arma::mat computeVolumeJacobian(double x, double y) = 0;
    virtual arma::mat computeBoundaryJacobian(double x, double y) = 0;
    
    virtual arma::mat computeBoundaryValue(double x, double y);
    
    arma::vec boundaryIntegral(int i, const MeshFn &u, int deg);
    arma::vec volumeIntegral(int i, int deg);
    
    // functor data
    // needed during assemble/jacobian routines
    arma::vec phi, psi, psi_x, psi_y;
    arma::mat UMinus, UPlus, U, UNeighbor;
    int m, iPlus, iMinus, iPhi, iPsi, neighbor;
    double nx, ny;
};
