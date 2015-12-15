#include <cmath>
#include <iostream>
#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
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
    
    PolyMesh msh;
    
    msh = quadUnitSquare(0.1);
    Advection eqn(msh);
    MeshFn f = MeshFn(msh, gaussian, 4);
    
    MeshFn b = eqn.assemble(f);
    
    f.gnuplot("plt/f.gnu");
    b.gnuplot("plt/b.gnu");
    
    return 0;
}
