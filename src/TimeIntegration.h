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
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct BackwardEuler : TimeStepper
{
    BackwardEuler(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct RK4 : TimeStepper
{
    RK4(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct RK2 : TimeStepper
{
    RK2(const MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { };
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};
