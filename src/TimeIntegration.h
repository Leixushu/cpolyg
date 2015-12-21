#pragma once

#include "MeshFn.h"
#include "MassMatrix.h"
#include "Equation.h"
#include "Timer/CH_Timer.H"

struct TimeStepper
{
    Equation &eqn;
    const MassMatrix &M;
    
    TimeStepper(const MassMatrix &a_M, Equation &a_eqn)
    : eqn(a_eqn), M(a_M) { };
    
    virtual MeshFn advance(const MeshFn &u, const double dt, const double t) = 0;
    virtual ~TimeStepper() {};
};

struct ForwardEuler : TimeStepper
{
    ForwardEuler(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0)
    {
        return u + dt*M.solve(eqn.assemble(u, t));
    }
};

struct RK4 : TimeStepper
{
    RK4(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0)
    {
        CH_TIMERS("RK4Advance");
        MeshFn k1 = dt*M.solve(eqn.assemble(u, t));
        MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1, t + 0.5*dt));
        MeshFn k3 = dt*M.solve(eqn.assemble(u + 0.5*k2, t + 0.5*dt));
        MeshFn k4 = dt*M.solve(eqn.assemble(u + k3, t + dt));
        
        return u + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    }
};

struct RK2 : TimeStepper
{
    RK2(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0)
    {
        MeshFn k1 = dt*M.solve(eqn.assemble(u, t));
        MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1, t + 0.5*dt));
        
        return u + 0.5*(k1 + k2);
    }
};
