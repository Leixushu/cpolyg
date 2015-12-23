#include "TimeIntegration.h"

MeshFn ForwardEuler::advance(const MeshFn &u, const double dt, const double t)
{
    return u + dt*M.solve(eqn.assemble(u, t));
}

MeshFn RK4::advance(const MeshFn &u, const double dt, const double t)
{
    CH_TIMERS("RK4Advance");
    MeshFn k1 = dt*M.solve(eqn.assemble(u, t));
    MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1, t + 0.5*dt));
    MeshFn k3 = dt*M.solve(eqn.assemble(u + 0.5*k2, t + 0.5*dt));
    MeshFn k4 = dt*M.solve(eqn.assemble(u + k3, t + dt));
    
    return u + (k1 + 2*k2 + 2*k3 + k4)/6.0;
}

MeshFn RK2::advance(const MeshFn &u, const double dt, const double t)
{
    MeshFn k1 = dt*M.solve(eqn.assemble(u, t));
    MeshFn k2 = dt*M.solve(eqn.assemble(u + 0.5*k1, t + 0.5*dt));
    
    return u + 0.5*(k1 + k2);
}
