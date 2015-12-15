#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"

double c5(double x, double y)
{
    return 5 + x + y;
}

int main(int argc, char ** argv)
{
    PolyMesh msh;
    
    msh = quadUnitSquare(0.99);
    
    MeshFn fn = MeshFn(msh, c5, 2);
    fn.a.print();
    
    return 0;
}