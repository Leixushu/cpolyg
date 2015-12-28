#include "MassMatrix.h"
#include "Legendre.h"

using namespace arma;
using namespace std;

MassMatrix::MassMatrix(PolyMesh &m, int d) : msh(m), deg(d)
{
    int i, j, k;
    int N;
    double integ;
    ProductFunctor prod(msh);
    
    basisSize = (deg+1)*(deg+2)/2;
    bl = basisSize;
    nb = msh.np;
    n_rows = nb;
    
    mat block(basisSize, basisSize);
    
    prod.phi = zeros<vec>(basisSize);
    prod.psi = zeros<vec>(basisSize);
    
    prod.m = deg+1;
    
    N = msh.np*basisSize;
    
    for (i = 0; i < msh.np; i++)
    {
        prod.i = i;
        
        for (j = 0; j < basisSize; j++)
        {
            prod.phi[j] = 1.0;
            for (k = 0; k < basisSize; k++)
            {
                prod.psi[k] = 1.0;
                
                integ = msh.polygonIntegral(prod, i);
                
                block(j, k) = integ;
                block(k, j) = integ;
                
                prod.psi[k] = 0.0;
            }
            prod.phi[j] = 0.0;
        }
        
        blocks.push_back(block);
        colIndices.push_back(i);
        rowBlock.push_back(i);
    }
    
    rowBlock.push_back(i);
}

MeshFn MassMatrix::solve(const MeshFn &fn) const
{
    int component;
    int i;
    MeshFn result(msh, deg, fn.nc);
    
    for (component = 0; component < fn.nc; component++)
    {
        for (i = 0; i < msh.np; i++)
        {
            vec b = fn.a.slice(i).col(component);
            vec x = arma::solve(blocks[i], b);
            
            result.a.slice(i).col(component) = x;
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
            vec b = blocks[i]*fn.a.slice(i).col(component);
            
            result.a.slice(i).col(component) = b;
        }
    }
    
    return result;
}
