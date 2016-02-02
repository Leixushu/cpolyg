#include "TimeIntegration.h"
#include "Preconditioners.h"

#define kNewtonMaxIterations 100
#define kNewtonTolerance 5.e-14

using namespace std;

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

MeshFn BackwardEuler::advance(const MeshFn &u, const double dt, const double t)
{
    MeshFn unp1 = u;
    MeshFn r(u.msh, u.deg, u.nc);
    MeshFn b(u.msh, u.deg, u.nc);
    int k;
    
    // Newton solve
    cout << "Beginning Newton solve" << endl;
    for (k = 0; k < kNewtonMaxIterations; k++)
    {
        unp1.gnuplot("plt/newton" + to_string(k) + ".gnu");
        
        b = eqn.assemble(unp1, t + dt);
        r = M.dot(unp1 - u) - dt*b;
        
        r.gnuplot("plt/residual" + to_string(k) + ".gnu");
        b.gnuplot("plt/b" + to_string(k) + ".gnu");
        
        cout << "    Iteration number " << k << ", residual norm = "
             << r.L2Norm().max() << endl;
        
        if (r.L2Norm().max() < kNewtonTolerance)
        {
            cout << "    Converged to tolerance" << endl;
            break;
        }
        
        Jacobian B = eqn.jacobian(unp1, t + dt);
        B *= -dt;
        B += M;
        //BlockJacobi pc(B);
        BlockILU0 pc(B);
        
        unp1 -= B.solve(r, pc);
    }
    
    return unp1;
}
