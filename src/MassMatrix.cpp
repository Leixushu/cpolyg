#include "MassMatrix.h"
#include "Legendre.h"
#include "blas/blas.h"
#include "Timer/CH_Timer.H"

using namespace arma;
using namespace std;

// Allocate and compute mass matrix
MassMatrix::MassMatrix(PolyMesh &m, int d) : msh(m), deg(d)
{
    int i, j, k;
    double integ;
    ProductFunctor prod(msh);
    
    // each block stores is (basis size) x (basis size)
    bl = (deg+1)*(deg+2)/2;
    
    // M is block diagonal, so number of blocks is number of rows
    nb = msh.np;
    n_rows = nb;
    
    // we store the blocks, and also the LU factorization
    blocks.clear();
    LU.clear();
    ipvt.resize(nb, Col<int>(bl));
    
    
    mat block(bl, bl);
    
    prod.phi = zeros<vec>(bl);
    prod.psi = zeros<vec>(bl);
    prod.m = deg+1;
    
    for (i = 0; i < msh.np; i++)
    {
        prod.i = i;
        
        for (j = 0; j < bl; j++)
        {
            prod.phi[j] = 1.0;
            for (k = 0; k <= j; k++)
            {
                prod.psi[k] = 1.0;
                
                // compute the integral of the product of basis functions
                integ = msh.polygonIntegral(prod, i);
                
                // the mass matrix is always symmetric
                block(j, k) = integ;
                block(k, j) = integ;
                
                prod.psi[k] = 0.0;
            }
            prod.phi[j] = 0.0;
        }
        
        // add the new block on the diagonal
        blocks.push_back(block);
        colIndices.push_back(i);
        rowBlock.push_back(i);
        
        // LAPACK LU factorization
        LU.push_back(block);
        cdgetrf(bl, bl, LU[i].memptr(), bl, ipvt[i].memptr());
        
    }
    
    rowBlock.push_back(i);
}

MeshFn MassMatrix::solve(const MeshFn &fn)
{
    int component;
    int i;
    MeshFn result(msh, deg, fn.nc);
    
    CH_TIMERS("Mass matrix solve");
    
    // the mass matrix is identical for each component
    for (component = 0; component < fn.nc; component++)
    {
        // since the mass matrix is block diagonal, we need to invert each block
        for (i = 0; i < msh.np; i++)
        {
            result.a.slice(i).col(component) = fn.a.slice(i).col(component);
            // use the precomputed LU factorization and pivots
            cdgetrs('N', bl, 1, LU[i].memptr(), bl, ipvt[i].memptr(), 
                    result.a.memptr() + i*bl*fn.nc + component*bl, bl);
        }
    }
    
    return result;
}

MeshFn MassMatrix::dot(const MeshFn &fn) const
{
    int component;
    int i;
    MeshFn result(msh, deg, fn.nc);
    
    for (component = 0; component < fn.nc; component++)
    {
        for (i = 0; i < msh.np; i++)
        {
            result.a.slice(i).col(component) = blocks[i]*fn.a.slice(i).col(component);
        }
    }
    
    return result;
}
