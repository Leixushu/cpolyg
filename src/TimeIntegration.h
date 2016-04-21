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

void getDIRK(int nStages, arma::mat &A, arma::vec &b, arma::vec &c);
struct DIRK : TimeStepper
{
    arma::mat A;
    arma::vec b, c;
    int nStages;
    DIRK(int a_nStages, MassMatrix &a_M, Equation &a_eqn)
        : TimeStepper(a_M, a_eqn), nStages(a_nStages) {
        getDIRK(nStages, A, b ,c);
    }
    
    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
};

void getIRK(int nStages, arma::mat &A, arma::mat &Ainv, arma::vec &b,
            arma::vec &c);
struct IRK : TimeStepper
{
    arma::mat A, Ainv;
    arma::vec b, c;
    int nStages;
    IRK(int a_nStages, MassMatrix &a_M, Equation &a_eqn)
        : TimeStepper(a_M, a_eqn), nStages(a_nStages) {
        getIRK(nStages, A, Ainv, b, c);
    }

    MeshFn advance(const MeshFn &u, const double dt, const double t = 0);
    MeshFn newAdvance(const MeshFn &u, const double dt, const double t = 0);
};
