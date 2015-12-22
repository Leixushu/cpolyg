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

int main(int argc, char ** argv)
{
    using namespace std;
    using namespace arma;
    
    int deg = 1;
    PolyMesh msh = quadUnitSquare(0.3);
    msh.gnuplot();
    
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    //BlockMatrix B = BlockMatrix::diag(mat(M.matrix), (deg+1)*(deg+2)/2);
    
    vec b(msh.np*M.basisSize); b.randu();
    
    NoPreconditioner pc;
    
    vec x1 = M.solve(b);
    vec x2 = zeros<vec>(msh.np*M.basisSize);
    
    int maxIt = 100;
    double tol = 1.e-12;
    int m = 20;
    //B.gmres(b, x2, m, tol, maxIt, pc);
    
//     vec x1 = M.matrix*b;
//     vec x2 = b;
//     vec x3 = B.matvec(b);
//     B.matvec(x2.memptr());
//     
//     cout << (x1 - x2).max() << endl;
//     cout << (x1 - x3).max() << endl;
    
    //cout << x1 - x2 << endl;
    
    cout << x1 - x2 << endl << endl << endl;
    
    cout << "Tol:  " << tol << endl;
    cout << "Iter: " << maxIt << endl;
    
    
    return 0;
}
