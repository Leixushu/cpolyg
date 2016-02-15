#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <cassert>
#include <iomanip>
#include "MeshFn.h"
#include "Quadrature.h"
#include "Legendre.h"
#include "Timer/CH_Timer.H"

using namespace std;
using namespace arma;

// Allocate a MeshFn with the given degree polynomial basis and given number of components
MeshFn::MeshFn(PolyMesh &a_msh, int a_deg, int a_nc /* = 1 */) : msh(a_msh)
{
    nc = a_nc;
    deg = a_deg;
    int basisSize = (deg+1)*(deg+2)/2;
    
    a = cube(basisSize, nc, msh.np);
}

// Create a (scalar) MeshFn by interpolating the given function with degree a_deg 
// polynomials
MeshFn::MeshFn(PolyMesh &a_msh, FnCallback cb, int a_deg) : msh(a_msh)
{
    FnCallbackFunctor functor(cb);
    int basisSize;
    nc = 1;
    deg = a_deg;
    
    basisSize = (deg+1)*(deg+2)/2;
    
    a = cube(basisSize, nc, msh.np);
    
    interp(functor, 0);
}

// copy constructor
MeshFn::MeshFn(const MeshFn &fn) : msh(fn.msh)
{
    deg = fn.deg;
    a = fn.a;
    nc = fn.nc;
}

// interpolate a scalar function given by a functor
void MeshFn::interp(const FnFunctor &cb, int component/* = 0 */)
{
    VecFnFunctor vecFunctor(cb);
    
    interp(vecFunctor, component);
}

// interpolate a vector function given by a functor
void MeshFn::interp(const VecFunctor &cb, int component/* = 0 */)
{
    gsl_integration_glfixed_table *glpts;
    vec quadraturePts(deg+1);
    int m;
    int basisSize;
    int p;
    int i, j, k;
    int numQuadPts;
    double ax,by,w,h;
    mat vals;
    mat G;
    
    // the function that we're interpolating should be scalar or vector -- not matrix
    // so we enforce one column only in the functor
    assert(cb.n_cols == 1);
    
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
    
    vals = mat(numQuadPts, cb.n_rows);
    
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
                vals.row(i*(deg+1) + j) = 
                    cb((quadraturePts[i]+1)*w/2.0 + ax,
                       (quadraturePts[j]+1)*h/2.0 + by).t();
            }
        }
        
        for (i = 0; i < cb.n_rows; i++)
        {
            a.slice(p).col(component + i) = solve(G, vals.col(i));
        }
    }
}

// Evaluate the given component of a MeshFn at x, y points in polynomial p
double MeshFn::eval(double x, double y, int p, int c/* = 0 */) const
{
    double xx, yy;
    vec coeffs;
    
    msh.getLocalCoordinates(p, x, y, xx, yy);
    
    coeffs = a.slice(p).col(c);
    
    return Leg2D(xx, yy, deg+1, coeffs);
}

// Write out a 'plot' of the function in a format readable by gnuplot 
void MeshFn::gnuplot(std::string filename) const
{
    int i, c;
    unsigned int j, k;
    int oldPrecision;
    double x1, x2, x3, y1, y2, y3, x, y, val;
    ofstream plotFile;
    mat &qr = Quadrature::tri2;
    
    CH_TIMERS("MeshFn gnuplot");
    
    plotFile.open(filename);
    
    oldPrecision = plotFile.precision();
    
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
            
            for (k = 0; k < qr.n_rows; k++)
            {
                plotFile.precision(oldPrecision);
                
                x = x1*(1-qr(k,0)-qr(k,1))+x2*qr(k,0)+x3*qr(k,1);
                y = y1*(1-qr(k,0)-qr(k,1))+y2*qr(k,0)+y3*qr(k,1);
                plotFile << x << "\t" << y;
                
                plotFile.precision(20);
                for (c = 0; c < nc; c++)
                {
                    val = eval(x, y, i, c);
                    plotFile << "\t" << val;
                }
                plotFile << endl;
            }
        }
    }
    
    plotFile.close();
}

MeshFn MeshFn::operator+(const MeshFn &fn2) const
{
    assert(deg == fn2.deg);
    assert(nc == fn2.nc);
    
    MeshFn fn(msh, deg, nc);
    fn.a = a + fn2.a;
    
    return fn; 
}

MeshFn MeshFn::operator-(const MeshFn &fn2) const
{
    assert(deg == fn2.deg);
    assert(nc == fn2.nc);
    
    MeshFn fn(msh, deg, nc);
    fn.a = a - fn2.a;
    
    return fn; 
}

MeshFn MeshFn::operator*(const double scale) const
{
    MeshFn result = *this;
    result.a *= scale;
    
    return result;
}

MeshFn & MeshFn::operator=(const MeshFn &fn)
{
    msh = fn.msh;
    deg = fn.deg;
    a = fn.a;
    nc = fn.nc;
    
    return *this;
}

mat MeshFn::L2Difference::operator()(double x, double y) const
{
    int c;
    vec error(fn.nc);
    vec exactValue = exact(x, y);
    
    for (c = 0; c < fn.nc; c++)
    {
        error(c) = pow(fn.eval(x, y, i, c) - exactValue(c), 2);
    }
    
    return error;
}

mat MeshFn::L2Functor::operator()(double x, double y) const
{
    int c;
    vec result(fn.nc);
    
    for (c = 0; c < fn.nc; c++)
    {
        result(c) = pow(fn.eval(x, y, i, c), 2);
    }
    
    return result;
}

double MeshFn::L2Error(const FnFunctor &exact) const
{
    VecFnFunctor vecExact(exact);
    
    return L2Error(vecExact)(0);
}

vec MeshFn::L2Error(const VecFunctor &exact) const
{
    int i;
    vec error = zeros<vec>(nc);
    L2Difference differenceSquared(*this, exact);
    
    for (i = 0; i < msh.np; i++)
    {
        differenceSquared.i = i;
        error += Quadrature::polygonIntegral(msh, differenceSquared, i, deg*2);
    }
    
    return sqrt(error);
}

// return the L2 norm of the function (as a vector of the L2 norm of the components)
vec MeshFn::L2Norm() const
{
    L2Functor functor(*this);
    vec l2 = zeros<vec>(nc);
    int i;
    
    for (i = 0; i < msh.np; i++)
    {
        functor.i = i;
        l2 += Quadrature::polygonIntegral(msh, functor, i, deg*2);
    }
    
    return sqrt(l2);
}
