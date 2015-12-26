#include "Meshes.h"
#include <vector>
#include <array>

using namespace std;

PolyMesh quadUnitSquare(double h)
{
    return quadRectangle(h, 1, 1);
}

PolyMesh hexUnitSquare(double h)
{
    return hexRectangle(h, 1, 1);
}

PolyMesh triUnitSquare(double h)
{
    return triRectangle(h, 1, 1);
}

PolyMesh honeycombUnitSquare(double h)
{
    return honeycombRectangle(h, 1, 1);
}

PolyMesh hexRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    int yi;

    for(x = 0; x <= width + kEPS; x += h)
    {
        yi = 0;
        for(y = 0; y <= height + kEPS; y += h)
        {
            if (x + 0.5*h*yi <= width + kEPS)
            {
                pt[0] = x + 0.5*h*yi;
                pt[1] = y;
                generatingPoints.push_back(pt);
            }
            
            yi = !yi;
        }
    }
    
    return PolyMesh(generatingPoints, width, height);
}


PolyMesh quadRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    
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
}

PolyMesh triRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
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
    
    msh.triangulate(generatingPoints);
    return msh;
}


PolyMesh honeycombRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    int yi;
    PolyMesh msh;

    for(x = 0; x <= width; x += h)
    {
        yi = 0;
        for(y = 0; y <= height; y += h)
        {
            if (x + 0.5*h*yi < width)
            {
                pt[0] = x + 0.5*h*yi;
                pt[1] = y;
                generatingPoints.push_back(pt);
            }
            
            yi = !yi;
        }
    }
    
    msh.triangulate(generatingPoints);
    return msh;
}
