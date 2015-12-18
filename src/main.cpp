#include <cmath>
#include <iostream>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include "MassMatrix.h"
#include "Advection.h"
#include "TimeIntegration.h"

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
    double h = 0.025;
    
    PolyMesh msh = hexUnitSquare(h);
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    Advection eqn(msh);
    
    MeshFn f = MeshFn(msh, gaussian, deg);
    f.gnuplot("plt/f.gnu");
    
    MeshFn unp1 = f;
    
    RK4 ti(M, eqn);
    
    int K;
    int i;
    double dt = h/30.0;
    
    K = M_PI/dt;
    K = 40;
    
    cout << "Computing total of " << K << " timesteps." << endl;
    
    for (i = 0; i < K; i++)
    {
        unp1 = ti.advance(unp1, dt);
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        
        if (i%10 == 0)
        {
            cout << "Timestep " << i << endl;
        }
    }
    
    return 0;
}
