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
    int yi;
    
    yi = 0;
    for(x = 0.5*h; x <= width + kEPS; x += h)
    {
        for(y = 0.5*h*yi; y <= height + kEPS; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            
            if (y > 0.5*h && y < height - h)
                pt[1] = y + h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            
            if (x > 0.5*h && x < width - h)
                pt[0] = x + h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            
            generatingPoints.push_back(pt);
        }
        if (x < width - 2*h) yi = !yi;
    }
    
    return PolyMesh(generatingPoints, width, height);
}

PolyMesh perturbedTriRectangle(double h, double p, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    int yi;
    
    ofstream ptfile;
    ptfile.open("plt/points.gnu");
    
    yi = 0;
    for(x = 0.5*h; x <= width + kEPS; x += h)
    {
        for(y = 0.5*h*yi; y <= height + kEPS; y += h)
        {
            pt[0] = x;
            pt[1] = y;
            
            if (y > 0.5*h && y < height - h)
                pt[1] = y + h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            
            if (x > 0.5*h && x < width - h)
                pt[0] = x + h*p*double(rand() - RAND_MAX/2)/RAND_MAX;
            
            generatingPoints.push_back(pt);
            ptfile << pt[0] << "\t\t" << pt[1] << endl;
        }
        if (x < width - 2*h) yi = !yi;
    }
    
    ptfile.close();
    
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

PolyMesh periodicRectangle(double h, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    PolyMesh msh;
    vector<int> polygon;
    int idx1, idx2, idx3, idx4;
    int Nx, Ny;
    int ix, iy, i, bci;
    PolyMesh::BoundaryInfo bc;
    
    bc.periodic = true;
    
    Nx = width/h;
    Ny = height/h;
    
    for(iy = 0; iy < Ny; iy++)
    {
        for(ix = 0; ix < Nx; ix++)
        {
        
            x = ix*h;
            y = iy*h;
            
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
            
            polygon = {idx1, idx2, idx3, idx4};
            msh.p.push_back(polygon);
        }
    }
    msh.np = msh.p.size();
    msh.p2p.resize(msh.np);
    
    i = 0;
    bci = -1;
    for(iy = 0; iy < Ny; iy++)
    {
        for(ix = 0; ix < Nx; ix++)
        {
            // bottom edge
            if (iy == 0)
            {
                bc.p1 = i;
                bc.p2 = Nx*(Ny-1) + ix;
                bc.a1 = msh.findVertex(ix*h, 0);
                bc.b1 = msh.findVertex((ix+1)*h, 0);
                bc.a2 = msh.findVertex(ix*h, Ny*h);
                bc.b2 = msh.findVertex((ix+1)*h, Ny*h);
                
                msh.bc.emplace(bci, bc);
                msh.p2p[i].push_back(bci);
                bci--;
            } else
            {
                msh.p2p[i].push_back(i - Nx);
            }
            
            // right edge
            if (ix == Nx-1)
            {
                bc.p1 = i;
                bc.p2 = i - Nx + 1;
                bc.a1 = msh.findVertex(Nx*h, iy*h);
                bc.b1 = msh.findVertex(Nx*h, (iy+1)*h);
                bc.a2 = msh.findVertex(0, iy*h);
                bc.b2 = msh.findVertex(0, (iy+1)*h);
                
                msh.bc.emplace(bci, bc);
                msh.p2p[i].push_back(bci);
                bci--;
            } else
            {
                msh.p2p[i].push_back(i+1);
            }
            
            // top edge
            if (iy == Ny-1)
            {
                bc.p1 = i;
                bc.p2 = ix;
                bc.a1 = msh.findVertex(ix*h, Ny*h);
                bc.b1 = msh.findVertex((ix+1)*h, Ny*h);
                bc.a2 = msh.findVertex(ix*h, 0);
                bc.b2 = msh.findVertex((ix+1)*h, 0);
                
                msh.bc.emplace(bci, bc);
                msh.p2p[i].push_back(bci);
                bci--;
            } else
            {
                msh.p2p[i].push_back(i+Nx);
            }
            
            // left edge
            if (ix == 0)
            {
                bc.p1 = i;
                bc.p2 = i + Nx - 1;
                bc.a1 = msh.findVertex(0, iy*h);
                bc.b1 = msh.findVertex(0, (iy+1)*h);
                bc.a2 = msh.findVertex(Nx*h, iy*h);
                bc.b2 = msh.findVertex(Nx*h, (iy+1)*h);
                
                msh.bc.emplace(bci, bc);
                msh.p2p[i].push_back(bci);
                bci--;
            } else
            {
                msh.p2p[i].push_back(i-1);
            }
            
            i++;
        }
    }
    
//     for (i = 0; i < msh.np; i++)
//     {
//         cout << "Polygon " << i << ", neighbors =   ";
//         unsigned int k;
//         for (k = 0; k < msh.p2p[i].size(); k++)
//         {
//             cout << msh.p2p[i][k] << ", ";
//         }
//         cout << endl;
//     }
    
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

    for (y = 0; y <= height + kEPS; y += sqrt(3)*h/2)
    {
        for (x = yi*0.5*h; x <= width + kEPS; x += h)
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
