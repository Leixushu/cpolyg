#include "Meshes.h"
#include <vector>
#include <array>

using namespace std;

PolyMesh quadUnitSquare(double h)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2>> generatingPoints;
    
    for(x = 0; x <= 1; x += h)
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
    
    for(x = 0; x <= 1; x += h)
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

PolyMesh triUnitSquare(double h)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2>> generatingPoints;
    PolyMesh msh;
    
    for(x = 0; x <= 1; x += h)
    {
        for(y = 0; y <= 1; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            generatingPoints.push_back(pt);
        }
    }
    
    msh.triangulate(generatingPoints);
    
    return msh;
}

PolyMesh quadRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2>> generatingPoints;
    PolyMesh msh;
    
    for(x = 0; x <= width; x += h)
    {
        for(y = 0; y <= height; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            generatingPoints.push_back(pt);
        }
    }
    
    return PolyMesh(generatingPoints, width, height);
    
    return msh;
}
