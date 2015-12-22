#include <cmath>
#include <iostream>
#include <iomanip>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include "MassMatrix.h"
#include "Advection.h"
#include "EulerVortex.h"
#include "TimeIntegration.h"
#include "BlockMatrix.h"
#include "Preconditioners.h"

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int deg = 1;
    PolyMesh msh = quadUnitSquare(0.3);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    BlockMatrix &B = M;
    
    vec b(msh.np*M.basisSize); b.randu();
    
    NoPreconditioner noPC;
    vec x1 = zeros<vec>(msh.np*M.basisSize);
    int maxIt = 100;
    double tol = 1.e-12;
    int m = 20;
    B.gmres(b, x1, m, tol, maxIt, noPC);
    cout << "Tol:  " << tol << endl;
    cout << "Iter: " << maxIt << endl;
    
    BlockJacobi PC(B);
    x1 = zeros<vec>(msh.np*M.basisSize);
    
    maxIt = 100;
    tol = 1.e-12;
    m = 20;
    B.gmres(b, x1, m, tol, maxIt, PC);
    
    cout << "Tol:  " << tol << endl;
    cout << "Iter: " << maxIt << endl;
    
    
    return 0;
}
