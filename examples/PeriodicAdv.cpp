#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "Advection.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

double gaussian(double x, double y)
{
    double x0, y0;
    
    x0 = 0.5 - 0.15;
    y0 = 0.5;
    
    return exp(-200*((x-x0)*(x-x0) + (y-y0)*(y-y0)));
}

int main(int argc, char ** argv)
{
    using namespace std;
    
    int deg = 2;
    double h = 0.1/2.0;
    
    PeriodicMesh msh = periodicRectangle(h, 1, 1);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    Advection eqn(msh);
    MeshFn unp1 = MeshFn(msh, gaussian, deg);
    
    int K = 2000;
    int i;
    double dt = 2.5/K;
    
    cout << "Using h = " << h << ", dt = " << dt << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    RK4 ti(M, eqn);
    
    for (i = 0; i < K; i++)
    {
        if (i%25 == 0)
        {
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
            unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        }
        
        unp1 = ti.advance(unp1, dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    return 0;
}
