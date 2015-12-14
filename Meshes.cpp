#include "Meshes.h"
#include <vector>
#include <array>

using namespace std;

PolyMesh oneSquare()
{
    PolyMesh msh;
    
    msh.np = 1;
    
    msh.v.resize(4);
    
    msh.v[0][0] = 0; msh.v[0][1] = 0;
    msh.v[1][0] = 1; msh.v[1][1] = 0;
    msh.v[2][0] = 1; msh.v[2][1] = 1;
    msh.v[3][0] = 0; msh.v[3][1] = 1;
    
    msh.p.resize(msh.np);
    msh.p[0].resize(4);
    msh.p[0][0] = 0;
    msh.p[0][1] = 1;
    msh.p[0][2] = 2;
    msh.p[0][3] = 3;
    
    msh.bb.resize(msh.np);
    msh.bb[0][0] = 0;
    msh.bb[0][1] = 0;
    msh.bb[0][2] = 1;
    msh.bb[0][3] = 1;
    
    return msh;
}

PolyMesh quadUnitSquare(double h)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2>> generatingPoints;
    
    for(x = 0; x <+ 1; x += h)
    {
        for(y = 0; y <= 1; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            generatingPoints.push_back(pt);
        }
    }
    
    return PolyMesh(generatingPoints);
}

PolyMesh hexUnitSquare(double h)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2>> generatingPoints;
    
    int yi;
    
    for(x = 0; x <+ 1; x += h)
    {
        yi = 0;
        for(y = 0; y <= 1; y += h)
        {
            pt[0] = x + 0.5*h*yi;
            pt[1] = y;
            generatingPoints.push_back(pt);
            
            yi = !yi;
        }
    }
    
    return PolyMesh(generatingPoints);
}