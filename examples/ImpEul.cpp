#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "EulerVortex.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

#define kGamma 1.4

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int i;
    int deg = 2;
    double h = 1.0;
    
    cout << "Using h = " << h << endl;
    
    PolyMesh msh = hexRectangle(h, 20, 15);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    EulerVortex eqn(msh, kGamma);
    
    MeshFn f = eqn.initialConditions(deg);
    MeshFn unp1 = f;
    
    f.gnuplot("plt/f.gnu");
    
    BackwardEuler ti(M, eqn);
    
    double dt = 1.0;
    
    int K = 100;
    
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
