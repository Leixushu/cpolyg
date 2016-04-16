#include <cmath>
#include <iostream>
#include <iomanip>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include "MassMatrix.h"
#include "EulerSetups.h"
#include "TimeIntegration.h"

#define kGamma 1.4

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int deg = 1;
    double h = 0.5;
    
    cout << "Using h = " << h << endl;
    
    //PolyMesh msh = hexRectangle(h, 20, 15);
    //PolyMesh msh = quadRectangle(h, 20, 15);
    //PolyMesh msh = periodicRectangle(100, 20, 20, 15);
    PolyMesh msh = periodicRectangle(100, 20, 20, 15);
    msh.gnuplot();
    
    VortexSolution exactSolution(kGamma);
    //BoundaryConditions bc = BoundaryConditions::dirichletConditions(msh, 
    //    &exactSolution);
    BoundaryConditions bc = BoundaryConditions::periodicConditions(msh);
    Euler eqn(msh, bc, kGamma);
    
    MassMatrix M(msh, deg, eqn.nc);
    M.spy("plt/M.gnu");
    
    MeshFn f = exactSolution.interpolated(msh, deg);
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
    K = 20.0/dt;
    
    cout << "Computing total of " << K << " timesteps." << endl;
    
    for (i = 0; i < K; i++)
    {
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        cout << "Beginning timestep " << i + 1 << endl;
        unp1 = ti.advance(unp1, dt, i*dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    exactSolution.t = i*dt;
    vec error = unp1.L2Error(exactSolution);
    cout << setprecision(20) << "L^2 errors: " << endl;
    cout << "Density:    " << error(0) << endl;
    cout << "Velocity u: " << error(1) << endl;
    cout << "Velocity v: " << error(2) << endl;
    cout << "Energy:     " << error(3) << endl;
    
    return 0;
}
