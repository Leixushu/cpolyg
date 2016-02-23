#include <cmath>
#include <iostream>
#include <iomanip>
#include "Meshes.h"
#include "EulerSetups.h"
#include "Jacobian.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

using namespace std;
using namespace arma;

#define kGamma 1.4

void solveit(int meshType, double cfl, int deg)
{
    int i;
    double h = 1.0;
    
    cout << "Creating mesh with h = " << h << endl;
    
    PolyMesh msh;
    
    cout << endl;
    switch(meshType)
    {
        case 0:
            cout << "### Hexagons ##############################################" << endl << endl;
            msh = hexRectangle(h/sqrt(6.0), 20, 15);
            break;
        case 1:
            cout << "### Squares ###############################################" << endl << endl;
            msh = quadRectangle(h*0.5*sqrt(sqrt(3)), 20, 15);
            break;
        case 2:
            cout << "### Right triangles #######################################" << endl << endl;
            msh = triRectangle(h*sqrt(sqrt(3)/2), 20, 15);
            break;
        case 3:
            cout << "### Equilateral triangles #################################" << endl << endl;
            msh = honeycombRectangle(h, 20, 15);
            break;
        case 4:
            cout << "### Peturbed polygons #####################################" << endl << endl;
            msh = perturbedQuadRectangle(h, 0.5, 20, 15);
            break;
        case 5:
            cout << "### Peturbed triangles ####################################" << endl << endl;
            msh = perturbedTriRectangle(h, 0.5, 20, 15);
            break;
    }
    
    cout << "Total degrees of freedom: " << msh.np*(deg+1)*(deg+2)*4/2 << endl;
    
    msh.gnuplot();
    
    VortexSolution exactSolution(kGamma);
    Euler eqn(msh, 
        BoundaryConditions::dirichletConditions(msh, &exactSolution), kGamma);
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    MeshFn unp1 = exactSolution.interpolated(msh, deg);
    
    BackwardEuler ti(M, eqn);
    
    double dt = h*cfl;
    int K = 1;
    cout << "Computing total of " << K << " timesteps." << endl;
    
    for (i = 0; i < K; i++)
    {
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
        cout << "Beginning timestep " << i + 1 << endl;
        unp1 = ti.advance(unp1, dt, i*dt);
    }
    unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
}

void doMeshes(double cfl, int deg)
{
    int i;
    
    cout << endl << endl << "CFL = " << cfl << endl << endl;
    
    for (i = 0; i < 4; i++)
    {
        solveit(i, cfl, deg);
    }
}

int main()
{
    doMeshes(1, 0);
    //doMeshes(1, 3);
    //doMeshes(2, 3);
    //doMeshes(4, 3);
}
