#include "MassMatrix.h"
#include "Legendre.h"
#include "Quadrature.h"
#include "blas/blas.h"
#include "Timer/CH_Timer.H"

using namespace arma;
using namespace std;

// Use the product functor to integrate the product of basis functions
struct ProductFunctor : FnFunctor
{
    arma::vec phi, psi;
    int i;
    int m;
    PolyMesh &msh;

    ProductFunctor(PolyMesh &a_msh) : msh(a_msh) {};

    double operator()(double x, double y) const
    {
        double xx, yy;
        msh.getLocalCoordinates(i, x, y, xx, yy);
    
        return Leg2D(xx, yy, m, phi)*Leg2D(xx, yy, m, psi);
    };
};

// Allocate and compute mass matrix
MassMatrix::MassMatrix(PolyMesh &a_msh, int a_deg, int a_nc)
    : msh(a_msh), deg(a_deg), nc(a_nc)
{
    int i, j, k;
    double integ;
    ProductFunctor prod(msh);
    
    int basisSize = (deg+1)*(deg+2)/2;
    
    // each block stores is (basis size) x (basis size)
    bl = nc*(deg+1)*(deg+2)/2;
    
    // M is block diagonal, so number of blocks is number of rows
    nb = msh.np;
    n_rows = nb;
    
    // we store the blocks, and also the LU factorization
    blocks.clear();
    LU.clear();
    ipvt.resize(nb, Col<int>(bl));
    
    
    mat block(bl, bl);
    
    prod.phi = zeros<vec>(basisSize);
    prod.psi = zeros<vec>(basisSize);
    prod.m = deg+1;
    
    for (i = 0; i < msh.np; i++)
    {
        prod.i = i;
        
        for (j = 0; j < basisSize; j++)
        {
            prod.phi[j] = 1.0;
            for (k = 0; k <= j; k++)
            {
                prod.psi[k] = 1.0;
                
                // compute the integral of the product of basis functions
                //integ = msh.polygonIntegral(prod, i);
                integ = Quadrature::polygonIntegral(msh, prod, i, deg*2);
                
                // the mass matrix is always symmetric
                for (int c = 0; c < nc; c++) {
                    block(j + c*basisSize, k + c*basisSize) = integ;
                    block(k + c*basisSize, j + c*basisSize) = integ;
                }
                
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
    int i;
    MeshFn result(msh, deg, fn.nc);
    
    CH_TIMERS("Mass matrix solve");
    
    // since the mass matrix is block diagonal, we need to invert each block
    for (i = 0; i < msh.np; i++)
    {
        result.a.slice(i) = fn.a.slice(i);
        // use the precomputed LU factorization and pivots
        cdgetrs('N', bl, 1, LU[i].memptr(), bl, ipvt[i].memptr(), 
                result.a.memptr() + i*bl*fn.nc, bl);
    }
    
    return result;
}

MeshFn MassMatrix::matvec(const MeshFn &fn) const
{
//     int component;
//     int i;
//     MeshFn result(msh, deg, fn.nc);
//     
//     for (component = 0; component < fn.nc; component++)
//     {
//         for (i = 0; i < msh.np; i++)
//         {
//             result.a.slice(i).col(component) = blocks[i]*fn.a.slice(i).col(component);
//         }
//     }
//     
//     return result;
    CH_TIMERS("Mass matrix matvec");
    MeshFn result = fn;
    BlockMatrix::matvec(result.a.memptr());
    
    return result;
}
