#include <cmath>
#include <iostream>
#include <iomanip>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include "MassMatrix.h"
#include "KelvinHelmholtz.h"
#include "TimeIntegration.h"

#define kGamma 1.4

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int deg = 2;
    double h = 0.025;
    
    cout << "Using h = " << h << endl;
    
    PolyMesh msh = periodicRectangle(h, 1, 1);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    KelvinHelmholtz eqn(msh, kGamma);
    
    MeshFn f = eqn.initialConditions(deg);
    MeshFn unp1 = f;
    
    //RK4 ti(M, eqn);
    RK2 ti(M, eqn);
    
    int K;
    int i;
    double dt = h/10.0;
    
    K = 1.0/dt;
    
    cout << "Computing total of " << K << " timesteps." << endl;
    for (i = 0; i < K; i++)
    {
        if (i%2 == 0)
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
