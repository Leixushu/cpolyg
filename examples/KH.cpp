#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "EulerSetups.h"
#include "TimeIntegration.h"

#define kGamma 1.4

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int deg = 2;
    double h = 0.01;
    
    cout << "Using h = " << h << endl;
    
    PolyMesh msh = periodicRectangle(400, 20, 1, 1);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    KelvinHelmholtz KH(kGamma);
    BoundaryConditions bc = BoundaryConditions::periodicConditions(msh);
    Euler eqn(msh, bc, kGamma);
    
    MeshFn unp1 = KH.interpolated(msh, deg);
    
    RK4 ti(M, eqn);
    //RK2 ti(M, eqn);
    
    int K;
    int i;
    double dt = h/20.0;
    cout << "Using dt = " << dt << endl;
    
    K = 1.0/dt;
    
    cout << "Computing total of " << K << " timesteps." << endl;
    for (i = 0; i < K; i++)
    {
        //if (i%2 == 0)
            unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        cout << "Beginning timestep " << i + 1 << endl;
        unp1 = ti.advance(unp1, dt, i*dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    return 0;
}
