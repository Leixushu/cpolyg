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
    
    int deg = 0;
    double h = 0.02;
    
    PolyMesh msh = honeycombUnitSquare(h);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    Advection eqn(msh);
    
    ExactGaussian exact(0);
    MeshFn f = MeshFn(msh, deg, 1);
    f.interp(exact);
    
    Jacobian B = eqn.jacobian(f);
    
    MeshFn unp1 = f;
    
    RK4 ti(M, eqn);
    
    int K;
    int i;
    double dt;
    
    K = 10*M_PI/h;
    //dt = M_PI/K/20;
    dt = 0.1;
    
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
