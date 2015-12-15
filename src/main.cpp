#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"
#include <cmath>

double c5(double x, double y)
{
    return 5 + x + y;
}

double gaussian(double x, double y)
{
    return exp(-150*((x-0.35)*(x-0.35) + (y-0.5)*(y-0.5)));
}

int main(int argc, char ** argv)
{
    PolyMesh msh;
    
    msh = quadUnitSquare(0.03);
    
    MeshFn fn = MeshFn(msh, gaussian, 2);
    
    fn.gnuplot("plt/test.gnu");
    
    return 0;
}
