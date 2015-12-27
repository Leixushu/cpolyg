#pragma once

#include <cstdlib>
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
    PolyMesh &msh;
    
    Equation(PolyMesh &m) : msh(m) { };
    
    virtual MeshFn assemble(const MeshFn &f, double t) = 0;
    virtual Jacobian jacobian(const MeshFn &f, double t)
    {
        std::cout << "Jacobian matrix not implemented for this equation" << std::endl;
        abort();
        return Jacobian(msh);
    };
    
    virtual ~Equation() {};
};

/**@}*/
