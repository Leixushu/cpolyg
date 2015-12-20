#include <cmath>
#include <iostream>
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
    
    int deg = 1;
    double h = 1;
    
    PolyMesh msh = quadRectangle(1, 20, 15);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    EulerVortex eqn(msh, kGamma);
    
    MeshFn f = eqn.initialConditions(deg);
    
    f.gnuplot("plt/u0.gnu");
    
//     MeshFn unp1 = f;
//     
//     RK4 ti(M, eqn);
//     
//     int K;
//     int i;
//     double dt = h/30.0;
//     
//     K = M_PI/dt;
//     
//     cout << "Computing total of " << K << " timesteps." << endl;
//     
//     for (i = 0; i < K; i++)
//     {
//         unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
//         unp1 = ti.advance(unp1, dt);
//         
//         if (i%4 == 0)
//         {
//             cout << "Timestep " << i << endl;
//         }
//     }
    
    return 0;
}
