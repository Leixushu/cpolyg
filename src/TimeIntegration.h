#pragma once

#include "MeshFn.h"
#include "MassMatrix.h"
#include "Equation.h"

struct TimeStepper
{
    Equation &eqn;
    const MassMatrix &M;
    
    TimeStepper(const MassMatrix &a_M, Equation &a_eqn)
    : eqn(a_eqn), M(a_M) { };
    virtual MeshFn advance(const MeshFn &u, const double dt) = 0;
};

struct ForwardEuler : TimeStepper
{
    ForwardEuler(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt)
    {
        return u + dt*M.solve(eqn.assemble(u));
    }
};

struct RK4 : TimeStepper
{
    RK4(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt)
    {
        MeshFn k1 = dt*M.solve(eqn.assemble(u));
        MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1));
        MeshFn k3 = dt*M.solve(eqn.assemble(u + 0.5*k2));
        MeshFn k4 = dt*M.solve(eqn.assemble(u + k3));
        
        return u + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    }
};
