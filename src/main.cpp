#include <cmath>
#include <iostream>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include "MassMatrix.h"
#include "Advection.h"

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
    
    int deg = 3;
    PolyMesh msh = quadUnitSquare(0.05);
    MassMatrix M(msh, deg);
    M.spy("plt/M.gnu");
    
    Advection eqn(msh);
    
    MeshFn f = MeshFn(msh, gaussian, deg);
    f.gnuplot("plt/f.gnu");
    
    MeshFn unp1 = f;
    
    int K;
    int i;
    double dt = 0.03;
    
    K = 600;
    
    for (i = 0; i < K; i++)
    {
        unp1 = unp1 + dt*eqn.assemble(unp1);
        unp1.gnuplot("plt/u" + to_string(i) + ".gnu");
    }
    
    return 0;
}
