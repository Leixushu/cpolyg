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

PolyMesh rectangle1D(double h, double width, double height) {
    double x;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    
    for(x = 0.5*h; x <= width + kEPS; x += h)
    {
        pt[0] = x;
        pt[1] = height/2;
        generatingPoints.push_back(pt);
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

PolyMesh periodicRectangle(int Nx, int Ny, double width, double height)
{
    double x, y;
    array<double, 2> pt;
    vector<array<double, 2> > generatingPoints;
    PolyMesh msh;
    vector<int> polygon;
    int idx1, idx2, idx3, idx4;
    int ix, iy, i, bci;
    double hx, hy;
    PolyMesh::BoundaryInfo bc;
    
    hx = width/Nx;
    hy = width/Ny;
    
    for(iy = 0; iy < Ny; iy++)
    {
        for(ix = 0; ix < Nx; ix++)
        {
        
            x = ix*hx;
            y = iy*hy;
            
            pt[0] = x;
            pt[1] = y;
            idx1 = msh.addVertex(pt);
            
            pt[0] = x+hx;
            pt[1] = y;
            idx2 = msh.addVertex(pt);
            
            pt[0] = x+hx;
            pt[1] = y+hy;
            idx3 = msh.addVertex(pt);
            
            pt[0] = x;
            pt[1] = y+hy;
            idx4 = msh.addVertex(pt);
            
            polygon = {idx1, idx2, idx3, idx4};
            msh.p.push_back(polygon);
        }
    }
    msh.np = msh.p.size();
    msh.p2p.resize(msh.np);
    
    bc.periodic = true;
    
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
                bc.a1 = msh.findVertex(ix*hx, 0);
                bc.b1 = msh.findVertex((ix+1)*hx, 0);
                bc.a2 = msh.findVertex(ix*hx, Ny*hy);
                bc.b2 = msh.findVertex((ix+1)*hx, Ny*hy);
                
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
                bc.a1 = msh.findVertex(Nx*hx, iy*hy);
                bc.b1 = msh.findVertex(Nx*hx, (iy+1)*hy);
                bc.a2 = msh.findVertex(0, iy*hy);
                bc.b2 = msh.findVertex(0, (iy+1)*hy);
                
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
                bc.a1 = msh.findVertex(ix*hx, Ny*hy);
                bc.b1 = msh.findVertex((ix+1)*hx, Ny*hy);
                bc.a2 = msh.findVertex(ix*hx, 0);
                bc.b2 = msh.findVertex((ix+1)*hx, 0);
                
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
                bc.a1 = msh.findVertex(0, iy*hy);
                bc.b1 = msh.findVertex(0, (iy+1)*hy);
                bc.a2 = msh.findVertex(Nx*hx, iy*hy);
                bc.b2 = msh.findVertex(Nx*hx, (iy+1)*hy);
                
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


/******************************************************************************/
#pragma mark -
/******************************************************************************/

PolyMesh naturalQuadOrdered(double h, double width, double height) {
    double x, y;
    array<double, 2> pt;
    PolyMesh msh;
    vector<int> polygon;
    
    int idx1, idx2, idx3, idx4;
    
    for(y = 0; y + h <= height + kEPS; y += h) {
        for(x = 0; x + h <= width + kEPS; x += h) {
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
    msh.computep2p();
    msh.computebb();
    msh.computeTriangulation();
    
    return msh;
}

PolyMesh naturalTriOrdered(double h, double width, double height) {
    double x, y;
    array<double, 2> pt;
    PolyMesh msh;
    vector<int> polygon;
    
    int idx1, idx2, idx3, idx4;
    
    for(y = 0; y + h <= height + kEPS; y += h) {
        for(x = 0; x + h <= width + kEPS; x += h) {
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
            
            polygon = {idx1, idx3, idx4};
            msh.p.push_back(polygon);
            polygon = {idx1, idx2, idx3};
            msh.p.push_back(polygon);
        }
    }
    
    msh.np = msh.p.size();
    msh.computep2p();
    msh.computebb();
    msh.computeTriangulation();
    
    return msh;
}

PolyMesh naturalHoneycombOrdered(double h, double width, double height) {
    double x, y;
    array<double, 2> pt;
    PolyMesh msh;
    vector<int> polygon;
    double l = sqrt(3)*h/2;
    
    int idx1, idx2, idx3, idx4;
    
    double xshift = 0;
    
    for(y = 0; y + l <= height + kEPS; y += l) {
        for(x = 0; x + h <= width + kEPS; x += h) {
            pt[0] = x + xshift;
            pt[1] = y;
            idx1 = msh.addVertex(pt);
            
            pt[0] = x+h + xshift;
            pt[1] = y;
            idx2 = msh.addVertex(pt);
            
            pt[0] = x+h/2 + xshift;
            pt[1] = y+l;
            idx3 = msh.addVertex(pt);
            
            pt[0] = x-h/2 + xshift;
            pt[1] = y+l;
            idx4 = msh.addVertex(pt);
            
            polygon = {idx1, idx3, idx4};
            msh.p.push_back(polygon);
            polygon = {idx1, idx2, idx3};
            msh.p.push_back(polygon);
        }
        xshift = -h/2*(xshift == 0);
    }
    
    msh.np = msh.p.size();
    msh.computep2p();
    msh.computebb();
    msh.computeTriangulation();
    
    return msh;
}

PolyMesh naturalHexOrdered(double h, double width, double height) {
    double x, y;
    array<double, 2> pt;
    PolyMesh msh;
    vector<int> polygon;
    double l = sqrt(3)*h;
    
    int idx1, idx2, idx3, idx4, idx5, idx6;
    
    for(y = 0; y + l <= height + kEPS; y += l) {
        double yshift = 0;
        for(x = 0; x + h <= width + kEPS; x += 3*h/2) {
            pt[0] = x + h;
            pt[1] = y + yshift;
            idx1 = msh.addVertex(pt);
            
            pt[0] = x + h/2;
            pt[1] = y + sqrt(3)*h/2 + yshift;
            idx2 = msh.addVertex(pt);
            
            pt[0] = x - h/2;
            pt[1] = y + sqrt(3)*h/2 + yshift;
            idx3 = msh.addVertex(pt);
            
            pt[0] = x - h;
            pt[1] = y + yshift;
            idx4 = msh.addVertex(pt);
            
            pt[0] = x - h/2;
            pt[1] = y - sqrt(3)*h/2 + yshift;
            idx5 = msh.addVertex(pt);
            
            pt[0] = x + h/2;
            pt[1] = y - sqrt(3)*h/2 + yshift;
            idx6 = msh.addVertex(pt);
            
            polygon = {idx1, idx2, idx3, idx4, idx5, idx6};
            msh.p.push_back(polygon);
            
            yshift = -sqrt(3)*h/2*(yshift == 0);
        }
    }
    
    msh.np = msh.p.size();
    msh.computep2p();
    msh.computebb();
    msh.computeTriangulation();
    
    return msh;
}
