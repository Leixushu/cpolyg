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
        
        //x0 = 0.5 + theta;
        //y0 = 0.5;
        
        return exp(-150*((x-x0)*(x-x0) + (y-y0)*(y-y0)));
    }
    
    ExactGaussian(double a_theta) : theta(a_theta) {};
};

int main(int argc, char ** argv)
{
    using namespace std;
    
    int deg = 4;
    double h = 0.1/2.0;
    
    double h2 = h*sqrt(2.0/3.0/sqrt(3.0));
    PolyMesh msh = hexUnitSquare(h2);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    Advection eqn(msh);
    
    ExactGaussian exact(0);
    MeshFn f = MeshFn(msh, deg, 1);
    f.interp(exact);
    
    MeshFn unp1 = f;
    
    int K = 1000;
    int i;
    double dt = M_PI/K;
    
    cout << "Using h = " << h << ", dt = " << dt << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    Jacobian B = eqn.jacobian(f, 0);
    
    for (i = 0; i < K; i++)
    {
        if (i%25 == 0)
        {
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
            unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        }
        
        MeshFn k1 = dt*M.solve(B.dot(unp1));
        MeshFn k2 = dt*M.solve(B.dot(unp1 + 0.5*k1));
        MeshFn k3 = dt*M.solve(B.dot(unp1 + 0.5*k2));
        MeshFn k4 = dt*M.solve(B.dot(unp1 + k3));
        
        unp1 += (k1 + 2*k2 + 2*k3 + k4)/6;
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    cout << setprecision(20) << "Computed until final time t=" << i*dt << endl;
    
    double l2err;
    exact.theta = i*dt;
    
    MeshFn theExact = MeshFn(msh, deg, 1);
    theExact.interp(exact);
    theExact.gnuplot("plt/exact.gnu");
    
    l2err = unp1.L2Error(exact);
    cout << "L^2 error = " << l2err << endl;
    
    return 0;
}
