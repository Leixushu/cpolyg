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
    int N = 1/h;
    
    PolyMesh msh = periodicRectangle(N, N, 1, 1);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    BoundaryConditions bc = BoundaryConditions::periodicConditions(msh);
    
    ConstantAdvection eqn(msh, bc, -0.75, 0.25);
    MeshFn unp1 = MeshFn(msh, gaussian, deg);
    
    int K = 1000;
    int i;
    double dt = 2.5/K;
    
    cout << "Using h = " << h << ", dt = " << dt << endl;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    Jacobian B = eqn.jacobian(unp1, 0);
    B *= -dt;
    B += M;
    BlockILU0 pc(B);
    
//     RK4 ti(M, eqn);
//     
    for (i = 0; i < K; i++)
    {
        if (i%25 == 0)
        {
            cout << "Beginning timestep " << i << ", t=" << i*dt << endl;
            unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        }
        
        unp1 = B.solve(M.dot(unp1), pc, kGMRESSolver);
        //unp1 = ti.advance(unp1, dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    
    return 0;
}
