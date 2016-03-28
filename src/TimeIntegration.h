#pragma once

#include "MeshFn.h"
#include "MassMatrix.h"
#include "Equation.h"
#include "Timer/CH_Timer.H"

// Abstract class for time stepping
struct TimeStepper
{
    Equation &eqn;
    MassMatrix &M;
    
    TimeStepper(MassMatrix &a_M, Equation &a_eqn) : eqn(a_eqn), M(a_M) { }
    
    virtual MeshFn advance(const MeshFn &u, const double dt, const double t) = 0;
    virtual ~TimeStepper() {};
};

struct ForwardEuler : TimeStepper
{
    ForwardEuler(MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct BackwardEuler : TimeStepper
{
    BackwardEuler(MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct RK4 : TimeStepper
{
    RK4(MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct RK2 : TimeStepper
{
    RK2(MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct DIRK3 : TimeStepper
{
    static arma::mat A, b, c;
    DIRK3(MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

struct IRK3 : TimeStepper
{
    static arma::mat A, b, c;
    IRK3(MassMatrix &a_M, Equation &a_eqn) : TimeStepper(a_M, a_eqn) { }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};
