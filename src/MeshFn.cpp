#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <cassert>
#include "MeshFn.h"
#include "Quadrature.h"
#include "Legendre.h"

using namespace std;
using namespace arma;

MeshFn::MeshFn(PolyMesh &a_msh, int a_deg, int a_nc) : msh(a_msh)
{
    nc = a_nc;
    deg = a_deg;
    int basisSize = (deg+1)*(deg+2)/2;
    
    a = cube(msh.np, nc, basisSize);
}

MeshFn::MeshFn(PolyMesh &a_msh, FnCallback cb, int a_deg) : msh(a_msh)
{
    int basisSize;
    nc = 1;
    deg = a_deg;
    
    basisSize = (deg+1)*(deg+2)/2;
    
    a = cube(msh.np, nc, basisSize);
    
    interp(cb, 0);
}

MeshFn::MeshFn(const MeshFn &fn) : msh(fn.msh)
{
    deg = fn.deg;
    a = fn.a;
    nc = fn.nc;
}

void MeshFn::interp(FnCallback cb, int component)
{
    gsl_integration_glfixed_table *glpts;
    vec quadraturePts(deg+1);
    int m;
    int basisSize;
    int p;
    int i, j, k;
    int numQuadPts;
    double ax,by,w,h;
    vec vals;
    mat G;
    
    // get the size of the basis given the degree
    basisSize = (deg+1)*(deg+2)/2;
    
    vec c = zeros<vec>(basisSize);
    
    // compute the quadrature points to use for the least squares interpolation
    m = deg+1;
    glpts = gsl_integration_glfixed_table_alloc(m);
    
    if (m % 2 == 0)
    {
        for (i = 0; i < m/2; i++)
        {
            quadraturePts[m/2 + i] = glpts->x[i];
            quadraturePts[m/2 - i - 1] = -glpts->x[i];
        }
    }
    else
    {
        quadraturePts[m/2] = 0.0;
        for (i = 1; i <= m/2; i++)
        {
            quadraturePts[m/2 - i] = -glpts->x[i];
            quadraturePts[m/2 + i] = glpts->x[i];
        }
    }
    // clean up
    gsl_integration_glfixed_table_free(glpts);
    
    // number of points in tensor product
    numQuadPts = (deg+1)*(deg+1);
    
    // create the matrix we use for the least squares
    G = mat(numQuadPts, basisSize);
    for (k = 0; k < basisSize; k++)
    {
        c[k] = 1.0;
        for (i = 0; i < deg + 1; i++)
        {
            for (j = 0; j < deg + 1; j++)
            {
                G(i*(deg+1) + j, k) = 
                      Leg2D(quadraturePts[i], quadraturePts[j], deg + 1, c);
            }
        }
        c[k] = 0.0;
    }
    
    vals = vec(numQuadPts);
    
    for (p = 0; p < msh.np; p++)
    {
        ax = msh.bb[p][0];
        by = msh.bb[p][1];
        w = msh.bb[p][2];
        h = msh.bb[p][3];
        
        for (i = 0; i < deg + 1; i++)
        {
            for (j = 0; j < deg + 1; j++)
            {
                vals(i*(deg+1) + j) = 
                    cb((quadraturePts[i]+1)*w/2.0 + ax, (quadraturePts[j]+1)*h/2.0 + by);
            }
        }
        
        a.tube(p, component) = solve(G, vals);
    }
}

double MeshFn::eval(double x, double y, int p, int c/* = 0 */)
{
    double xx, yy;
    vec coeffs;
    
    msh.getLocalCoordinates(p, x, y, xx, yy);
    
    //xx = 2.0*(x - msh.bb[p][0])/msh.bb[p][2] - 1;
    //yy = 2.0*(y - msh.bb[p][1])/msh.bb[p][3] - 1;
    
    coeffs = a.tube(p, c);
    
    return Leg2D(xx, yy, deg+1, coeffs);
}

void MeshFn::gnuplot(std::string filename, int c/* = 0 */)
{
    int i, j, k;
    double x1, x2, x3, y1, y2, y3, x, y, val;
    ofstream plotFile;
    auto &qr = Quadratures::tri2;
    
    plotFile.open(filename);
    
    for (i = 0; i < msh.np; i++)
    {
        for (j = 0; j < msh.tri[i].triangles.size(); j++)
        {
            x1 = msh.tri[i].points[msh.tri[i].triangles[j][0]][0];
            x2 = msh.tri[i].points[msh.tri[i].triangles[j][1]][0];
            x3 = msh.tri[i].points[msh.tri[i].triangles[j][2]][0];
            y1 = msh.tri[i].points[msh.tri[i].triangles[j][0]][1];
            y2 = msh.tri[i].points[msh.tri[i].triangles[j][1]][1];
            y3 = msh.tri[i].points[msh.tri[i].triangles[j][2]][1];
            
            for (k = 0; k < qr.size(); k++)
            {
                x = x1*(1-qr[k][0]-qr[k][1])+x2*qr[k][0]+x3*qr[k][1];
                y = y1*(1-qr[k][0]-qr[k][1])+y2*qr[k][0]+y3*qr[k][1];
                val = eval(x, y, i, c);
                plotFile << x << "\t" << y << "\t" << val << endl;
            }
        }
    }
    
    plotFile.close();
}

MeshFn MeshFn::operator+(MeshFn fn2)
{
    assert(deg == fn2.deg);
    assert(nc == fn2.nc);
    
    MeshFn fn(msh, deg, nc);
    fn.a = a + fn2.a;
    
    return fn; 
}

MeshFn & MeshFn::operator=(const MeshFn fn)
{
    msh = fn.msh;
    deg = fn.deg;
    a = fn.a;
    nc = fn.nc;
    
    return *this;
}

MeshFn operator*(const double scale, const MeshFn& fn)
{
    MeshFn result = fn;
    result.a *= scale;
    return result;
}
