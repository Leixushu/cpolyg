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

double planarWave(double x, double y)
{
    int nx = 10;
    int ny = 10;
    return sin(nx*x + ny*y);
}

int main(int argc, char ** argv)
{
    using namespace std;
    
    int deg = 0;
    double h = 0.02;
    
    PolyMesh msh = triUnitSquare(h);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    Advection eqn(msh);
    
    //FnCallbackFunctor exact(planarWave);
    
    ExactGaussian exact(0);
    
    MeshFn f = MeshFn(msh, deg, 1);
    f.interp(exact);
    
    Jacobian B = eqn.jacobian(f, 0);
    B.spy("plt/J.gnu");
    
    MeshFn unp1 = f;
    
    int K;
    int i;
    double dt;
    
    dt = 0.05;
    K = 100;
    
    cout << "Using h = " << h << ", dt = " << dt << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    B *= -dt;
    B += M;
    
    B.spy("plt/B.gnu");
    
    BlockJacobi pc(B);
    
    for (i = 0; i < K; i++)
    {
        if (i%10 == 0)
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
        
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        unp1 = B.solve(M.dot(unp1), pc, kJacobiSolver);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    cout << setprecision(20) << "Computed until final time t=" << i*dt << endl;
    
    double l2err;
    
    exact.theta = i*dt;
    l2err = unp1.L2Error(exact);
    cout << "L^2 error = " << l2err << endl;
    
    return 0;
}
