#include "PolyMesh.h"
#include "MeshFn.h"
#include "Meshes.h"

double c5(double x, double y)
{
    return 2*x - 1 + 2*y - 1;
}

int main(int argc, char ** argv)
{
    PolyMesh msh;
    
    msh = hexUnitSquare(0.1);
    
//     MeshFn fn = MeshFn(msh, c5, 2);
//     
//     fn.a.print();
    
    return 0;
}