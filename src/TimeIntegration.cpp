#include "TimeIntegration.h"
#include "Preconditioners.h"

#define kNewtonMaxIterations 100
#define kNewtonTolerance 5.e-10

using namespace std;
using namespace arma;

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
        //unp1.gnuplot("plt/newton" + to_string(k) + ".gnu");
        
        b = eqn.assemble(unp1, t + dt);
        r = M.matvec(unp1 - u) - dt*b;
        
        //r.gnuplot("plt/residual" + to_string(k) + ".gnu");
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
        BlockILU0 pc(B);
        
        unp1 -= B.solve(r, pc, kGMRESSolver);
    }
    
    return unp1;
}

static double const alpha=0.435866521508459;
mat DIRK3::A = {{(1+alpha)/2, 0,           0},
                {-alpha/2,    (1+alpha)/2, 0},
                {1+alpha,    -(1+2*alpha), (1+alpha)/2}};
mat DIRK3::b = {1.0/6.0/(alpha*alpha),1-1.0/3.0/(alpha*alpha), 1.0/6.0/(alpha*alpha)};
mat DIRK3::c = {(1+alpha)/2,1.0/2.0,(1-alpha)/2};

MeshFn DIRK3::advance(const MeshFn &u, const double dt, const double t)
{
    const static int nStages = 3;
    int stage;
    MeshFn zero(u.msh, u.deg, u.nc);
    zero.a.fill(0);
    array<MeshFn, nStages> k = {{zero, zero, zero}};
    
    for (stage = 0; stage < nStages; stage++)
    {
        cout << "DIRK stage " << stage+1 << endl;
        
        cout << "  Beginning Newton solve" << endl;
        for (int iter = 0; iter < kNewtonMaxIterations; iter++)
        {
            MeshFn cu = u;
            for (int j = 0; j <= stage; j++)
            {
                cu += A(stage, j)*k[j];
            }
            
            MeshFn rhs = eqn.assemble(cu, t + dt*c(stage));
            MeshFn r = M.matvec(k[stage]) - dt*rhs;
            
            cout << "    Iteration number " << iter << ", residual norm = "
                 << r.L2Norm().max() << endl;
            
            if (r.L2Norm().max() < kNewtonTolerance)
            {
                cout << "    Converged to tolerance" << endl;
                break;
            }
            
            Jacobian B = eqn.jacobian(cu, t + dt*c(stage));
            B *= -dt*A(stage, stage);
            B += M;
            BlockILU0 pc(B);
            k[stage] -= B.solve(r, pc, kGMRESSolver);
        }
    }
    
    MeshFn unp1 = u;
    
    for (int i = 0; i < nStages; i++)
    {
        unp1 += b(i)*k[i];
    }
    
    return unp1;
}
