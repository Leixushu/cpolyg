#include <cmath>
#include "Meshes.h"
#include "Advection.h"
#include "TimeIntegration.h"

// Initial conditions
double gaussian(double x, double y)
{
    return exp(-150*((x-0.25)*(x-0.25) + (y-0.5)*(y-0.5)));
}

int main(int argc, char ** argv)
{
    // specify the degree of polynomials to use
    int deg = 2;
    // specify the size of each element in the mesh
    double h = 0.02;
    // specify the time step
    double dt = 0.1;
    
    // create the mesh, in this case hexagonal tessellation
    // of the unit square
    PolyMesh msh = hexUnitSquare(h);
    
    // specify the advection equation
    Advection eqn(msh);
    
    // compute the mass matrix
    MassMatrix M(msh, deg);
    
    // create a mesh function by interpolating the 
    // initial conditions defined by the function above
    MeshFn f = MeshFn(msh, gaussian, deg);
    
    // create the time integration object
    BackwardEuler ti(M, eqn);
    
    for (int i = 0; i < 10; i++)
    {
        f = ti.advance(f, dt, i*dt);
        f.gnuplot("plt/u" + std::to_string(i) + ".gnu");
    }
    
    return 0;
}













// #include <cmath>
// #include <iostream>
// #include <iomanip>
// #include "PolyMesh.h"
// #include "MeshFn.h"
// #include "Meshes.h"
// #include "MassMatrix.h"
// #include "Advection.h"
// #include "EulerVortex.h"
// #include "TimeIntegration.h"
// #include "BlockMatrix.h"
// #include "Preconditioners.h"
// 
// int main(int argc, char ** argv)
// {
//     using namespace std;
//     using namespace arma;
//     
//     int deg = 1;
//     PolyMesh msh = quadUnitSquare(0.3);
//     msh.gnuplot();
//     
//     MassMatrix M(msh, deg);
//     M.spy("plt/M.gnu");
//     
//     BlockMatrix &B = M;
//     
//     vec b(msh.np*M.basisSize); b.randu();
//     
//     NoPreconditioner noPC;
//     vec x1 = zeros<vec>(msh.np*M.basisSize);
//     int maxIt = 100;
//     double tol = 1.e-12;
//     int m = 20;
//     B.gmres(b, x1, m, tol, maxIt, noPC);
//     cout << "Tol:  " << tol << endl;
//     cout << "Iter: " << maxIt << endl;
//     
//     BlockJacobi PC(B);
//     x1 = zeros<vec>(msh.np*M.basisSize);
//     
//     maxIt = 100;
//     tol = 1.e-12;
//     m = 20;
//     B.gmres(b, x1, m, tol, maxIt, PC);
//     
//     cout << "Tol:  " << tol << endl;
//     cout << "Iter: " << maxIt << endl;
//     
//     
//     return 0;
// }
