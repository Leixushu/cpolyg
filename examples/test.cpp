#include <cmath>
#include "Meshes.h"
#include "Advection.h"
#include "TimeIntegration.h"
#include "Preconditioners.h"

// Initial conditions
double gaussian(double x, double y)
{
    return exp(-150*((x-0.25)*(x-0.25) + (y-0.5)*(y-0.5)));
}

double zero(double x, double y)
{
    return 0;
}

int main(int argc, char ** argv)
{
    using namespace std;
    // specify the degree of polynomials to use
    int deg = 2;
    // specify the size of each element in the mesh
    double h = 0.1;
    // specify the time step
    double dt = h/10;
    
    // create the mesh, in this case hexagonal tessellation
    // of the unit square
    PolyMesh msh = hexUnitSquare(h);
    
    // create a mesh function by interpolating the 
    // initial conditions defined by the function above
    MeshFn f = MeshFn(msh, gaussian, deg);
    
    // specify the advection equation
    // with zero Dirichlet boundary conditions
    VecFnCallbackFunctor zeroFunctor(zero);
    BoundaryConditions bc = 
        BoundaryConditions::dirichletConditions(msh, &zeroFunctor);
    Advection eqn(msh, bc);
    
    // compute the mass matrix
    MassMatrix M(msh, deg);
    
    // create the time integration object
    ForwardEuler ti(M, eqn);
    
    for (int i = 0; i < 200; i++)
    {
        cout << "Beginning timestep " << i << endl;
        f = ti.advance(f, dt, i*dt);
        f.gnuplot("plt/u" + std::to_string(i) + ".gnu");
    }
    
    return 0;
}




//     cout << "Computing the Jacobian" << endl;
//     Jacobian J = eqn.jacobian(f, 0);
//     
//     cout << "Extracting the block diagonal" << endl;
//     BlockMatrix DJ = BlockMatrix::diag(J);
//     BlockMatrix OJ = BlockMatrix::offDiag(J);
//     
//     J.spy("plt/J.gnu");
//     DJ.spy("plt/DJ.gnu");
//     OJ.spy("plt/OJ.gnu");
//     
//     BlockILU0 ilu(J);
//     
//     ilu.DD.spy("plt/DD.gnu");
//     ilu.DU.spy("plt/DU.gnu");
//     ilu.L.spy("plt/L.gnu");








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
