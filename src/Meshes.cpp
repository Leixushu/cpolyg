#include "Meshes.h"
#include <vector>
#include <array>

using namespace std;

PolyMesh hexRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    int yi;

    for(x = 0; x <= width + kEPS; x += 3*h)
    {
        yi = 0;
        for(y = 0; y <= height + kEPS; y += 0.5*sqrt(3)*h)
        {
            if (x + 1.5*h*yi <= width + kEPS)
            {
                pt[0] = x + 1.5*h*yi;
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
    
    for(x = 0.5*h; x <= width + kEPS; x += h)
    {
        for(y = 0.5*h; y <= height + kEPS; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            generatingPoints.push_back(pt);
        }
    }
    
    return PolyMesh(generatingPoints, width, height);
}

PolyMesh perturbedQuadRectangle(double h, double p, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    
    for(x = 0.5*h; x <= width + kEPS; x += h)
    {
        for(y = 0.5*h; y <= height + kEPS; y += h)
        {
            pt[0] = x + h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            pt[1] = y + h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            generatingPoints.push_back(pt);
        }
    }
    
    return PolyMesh(generatingPoints, width, height);
}

PolyMesh perturbedTriRectangle(double h, double p, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    
    for(x = 0.5*h; x <= width + kEPS; x += h)
    {
        for(y = 0.5*h; y <= height + kEPS; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            
            if (x > 0.5*h && x < width - h && y > 0.5*h && y < height - h)
            {
                pt[0] += h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
                pt[1] += h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            }
            
            generatingPoints.push_back(pt);
        }
    }
    
    return PolyMesh::triangulate(generatingPoints);
}

PolyMesh triRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    PolyMesh msh;
    vector<int> polygon;
    int idx1, idx2, idx3, idx4;
    
    for(x = 0; x + h <= width + kEPS; x += h)
    {
        for(y = 0; y + h <= height + kEPS; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            idx1 = msh.addVertex(pt);
            
            pt[0] = x+h;
            pt[1] = y;
            idx2 = msh.addVertex(pt);
            
            pt[0] = x+h;
            pt[1] = y+h;
            idx3 = msh.addVertex(pt);
            
            pt[0] = x;
            pt[1] = y+h;
            idx4 = msh.addVertex(pt);
            
            polygon = {idx1, idx2, idx3};
            msh.p.push_back(polygon);
            polygon = {idx1, idx3, idx4};
            msh.p.push_back(polygon);
        }
    }
    
    msh.np = msh.p.size();
    msh.computep2p();
    msh.computebb();
    msh.computeTriangulation();
    
    return msh;
}


PolyMesh honeycombRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    int yi;

    yi = 0;
    for(x = 0; x <= width + kEPS; x += sqrt(3)*h/2)
    {
        for(y = yi*0.5*h; y <= height + kEPS; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            generatingPoints.push_back(pt);
        }
        
        yi = !yi;
    }
    
    return PolyMesh::triangulate(generatingPoints);
}

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
