#pragma once

#include <cstdlib>
#include <armadillo>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Jacobian.h"

/** \defgroup Equations Equations
 *
 *  @brief Implementations of discontinuous Galerkin discretizations of various PDEs
 *  @{
 */

/// Abstract class for implementing a PDE
struct Equation
{
    struct VolumeTermFunctor : VecFunctor
    {
        int m;
        int i;
        
        const arma::vec *psi_x, *psi_y;
        arma::mat U;
        PolyMesh &msh;
        
        VolumeTermFunctor(PolyMesh &a_msh) : msh(a_msh) { };
        
        virtual arma::mat operator()(double x, double y) const = 0;
    };
    
    struct NumericalFluxFunctor : VecFunctor
    {
        double nx, ny;
        int m;
        int iMinus;
        int iPlus;
        
        const arma::vec *psi;
        arma::mat UMinus, UPlus;
        
        PolyMesh &msh;
        
        NumericalFluxFunctor(PolyMesh &a_msh) : msh(a_msh) { };
        
        virtual arma::mat operator()(double x, double y) const = 0;
    };
    
    struct VolumeTermJacobianFunctor : VecFunctor
    {
        int m;
        int i;
        
        const arma::vec *psi_x, *psi_y, *phi;
        arma::mat U;
        
        PolyMesh &msh;
        
        VolumeTermJacobianFunctor(PolyMesh &a_msh) : msh(a_msh) { };
        
        virtual arma::mat operator()(double x, double y) const = 0;
    };
    
    struct NumericalFluxJacobianFunctor : VecFunctor
    {
        double nx, ny;
        int m;
        int iPhi, iPsi, neighbor;
        
        const arma::vec *psi, *phi;
        arma::mat U, UNeighbor;
        PolyMesh &msh;
        
        NumericalFluxJacobianFunctor(PolyMesh &a_msh) : msh(a_msh) { };
        
        virtual arma::mat operator()(double x, double y) const = 0;
    };
    
    PolyMesh &msh;
    
    VolumeTermFunctor *volumeTerm;
    NumericalFluxFunctor *boundaryTerm;
    
    VolumeTermJacobianFunctor *volumeJacobian;
    NumericalFluxJacobianFunctor *boundaryDerivative;
    
    int nc;
    
    Equation(PolyMesh &a_msh);
    
    arma::vec boundaryIntegral(int i, const arma::vec &psi, const MeshFn &U);
    arma::vec volumeIntegral(int i, const arma::vec &psi_x, const arma::vec &psi_y);
    
    MeshFn assemble(const MeshFn &f, double t);
    Jacobian jacobian(const MeshFn &f, double t);
};

/**@}*/
