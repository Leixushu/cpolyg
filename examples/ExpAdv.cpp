#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "Advection.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

struct ExactGaussian : FnFunctor
{
    double theta;
    
    double operator()(double x, double y) const
    {
        double x0, y0;
        
        x0 = 0.5 - 0.15*cos(2*theta);
        y0 = 0.5 + 0.3*cos(theta)*sin(theta);
        
        return exp(-150*((x-x0)*(x-x0) + (y-y0)*(y-y0)));
    }
    
    ExactGaussian(double a_theta) : theta(a_theta) {};
};

double c5(double x, double y)
{
    return 5 + x + 7*x*y + 3*x*x*y;
}

double gaussian(double x, double y)
{
    return exp(-150*((x-0.35)*(x-0.35) + (y-0.5)*(y-0.5)));
}

int main(int argc, char ** argv)
{
    using namespace std;
    
    int deg = 3;
    double h = 0.03;
    
    PolyMesh msh = hexUnitSquare(h);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    Advection eqn(msh);
    
    ExactGaussian exact(0);
    MeshFn f = MeshFn(msh, deg, 1);
    f.interp(exact);
    
    MeshFn unp1 = f;
    
    //RK4 ti(M, eqn);
    //ForwardEuler ti(M, eqn);
    
    int K = 500;
    int i;
    double dt = 1.0/K;
    
    cout << "Using h = " << h << ", dt = " << dt << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    Jacobian B = eqn.jacobian(f, 0);
    BlockILU0 pc(B);
    
    for (i = 0; i < K; i++)
    {
        if (i%25 == 0)
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
            unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        
        MeshFn k1 = dt*M.solve(B.dot(unp1));
        MeshFn k2 = dt*M.solve(B.dot(unp1 + 0.5*k1));
        MeshFn k3 = dt*M.solve(B.dot(unp1 + 0.5*k2));
        MeshFn k4 = dt*M.solve(B.dot(unp1 + k3));
        
        unp1 += (k1 + 2*k2 + 2*k3 + k4)/6;
        
        //unp1 = ti.advance(unp1, dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    cout << setprecision(20) << "Computed until final time t=" << i*dt << endl;
    
    double l2err;
    exact.theta = i*dt;
    l2err = unp1.L2Error(exact);
    cout << "L^2 error = " << l2err << endl;
    
    return 0;
}
