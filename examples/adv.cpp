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
#include "Timer/CH_Timer.H"

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
    
    int deg = 1;
    double h = 0.1;
    
    CH_TIMERS("AdvectionMain");
    
    PolyMesh msh = quadUnitSquare(h);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    Advection eqn(msh);
    
    MeshFn f = MeshFn(msh, gaussian, deg);
    MeshFn unp1 = f;
    
    RK4 ti(M, eqn);
    
    int K;
    int i;
    double dt;
    
    K = 8*M_PI/h;
    dt = M_PI/K;
    
    K = 50;
    
    cout << "Using h = " << h << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    for (i = 0; i < K; i++)
    {
        if (i%10 == 0)
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
        
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        unp1 = ti.advance(unp1, dt);
    }
    
    cout << setprecision(20) << "Computed until final time t=" << i*dt << endl;
    
    double l2err;
    FnCallbackFunctor exact(gaussian);
    
    l2err = unp1.L2Error(exact);
    cout << "L^2 error = " << l2err << endl;
    
    return 0;
}
