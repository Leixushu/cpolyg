#include <cmath>
#include <iostream>
#include <iomanip>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include "MassMatrix.h"
#include "Advection.h"
#include "EulerVortex.h"
#include "TimeIntegration.h"

#define kGamma 1.4

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
    using namespace arma;
    
    int deg = 1;
    double h = 0.75;
    
    cout << "Using h = " << h << endl;
    
    PolyMesh msh = hexRectangle(h, 20, 15);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    EulerVortex eqn(msh, kGamma);
    
    MeshFn f = eqn.initialConditions(deg);
    MeshFn unp1 = f;
    
    f.gnuplot("plt/f.gnu");
    
    //RK4 ti(M, eqn);
    RK2 ti(M, eqn);
    //ForwardEuler ti(M, eqn);
    
    int K;
    int i;
    double dt = h/30.0;
    
    //K = M_PI/dt;
    //K = sqrt(125)/dt;
    K = 1.0/dt;
    
    cout << "Computing total of " << K << " timesteps." << endl;
    
    for (i = 0; i < K; i++)
    {
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        cout << "Beginning timestep " << i + 1 << endl;
        unp1 = ti.advance(unp1, dt, i*dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    eqn.exact->t = i*dt;
    vec error = unp1.L2Error(*eqn.exact);
    cout << setprecision(20) << "L^2 errors: " << endl;
    cout << "Density:    " << error(0) << endl;
    cout << "Velocity u: " << error(1) << endl;
    cout << "Velocity v: " << error(2) << endl;
    cout << "Energy:     " << error(3) << endl;
    
    return 0;
}
