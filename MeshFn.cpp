#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include "MeshFn.h"
#include "leg.h"

using namespace std;
using namespace arma;

MeshFn::MeshFn(PolyMesh &a_msh, int a_deg, int a_nc) : msh(a_msh)
{
    int basis_size = (deg+1)*(deg+2)/2;
    nc = a_nc;
    deg = a_deg;
    
    a = cube(msh.np, nc, basis_size);
}

MeshFn::MeshFn(PolyMesh &a_msh, FnCallback cb, int a_deg) : msh(a_msh)
{
    int basis_size;
    nc = 1;
    deg = a_deg;
    
    basis_size = (deg+1)*(deg+2)/2;
    
    a = cube(msh.np, nc, basis_size);
    
    interp(cb, 0);
}

void MeshFn::interp(FnCallback cb, int component)
{
    gsl_integration_glfixed_table *glpts;
    int basis_size;
    int p;
    int i, j, k;
    int numQuadPts;
    double ax,by,w,h;
    vector<double> c;
    vec vals;
    mat G;
    
    // get the size of the basis given the degree
    basis_size = (deg+1)*(deg+2)/2;
    
    c.resize(basis_size, 0.0);
    
    // compute the quadrature points to use for the least squares interpolation
    glpts = gsl_integration_glfixed_table_alloc(deg+1);
    
    // number of points in tensor product
    numQuadPts = (deg+1)*(deg+1);
    
    // create the matrix we use for the least squares
    G = mat(numQuadPts, basis_size);
    for (k = 0; k < basis_size; k++)
    {
        c[k] = 1.0;
        for (i = 0; i < deg + 1; i++)
        {
            for (j = 0; j < deg + 1; j++)
            {
                G(i*(deg+1) + j, k) = 
                    Leg2D(glpts->x[i], glpts->x[j], deg + 1, c.data());
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
                    cb((glpts->x[i]+1)*w/2.0 + ax, (glpts->x[j]+1)*h/2.0 + by);
            }
        }
        
        a.tube(p, component) = solve(G, vals);
    }
    
    gsl_integration_glfixed_table_free(glpts);

}

PolyFn MeshFn::getPolyFn(int p, int c)
{
    PolyFn fn;
    fn.a = a.tube(p, c);
    
    return fn;
}