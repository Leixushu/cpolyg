#include "Equation.h"

using namespace arma;
using namespace std;

Equation::Equation(PolyMesh &a_msh) : msh(a_msh)
{
    volumeTerm = NULL;
    boundaryTerm = NULL;
    volumeJacobian = NULL;
    boundaryDerivative = NULL;
}

vec Equation::boundaryIntegral(int i, const vec &psi, const MeshFn &U)
{
    int a1, b1, a2, b2;
    int nv1, nv2;
    int j, l, i2;
    unsigned int k;
    int neighbor;
    vec integ = zeros<vec>(nc);
    mat UNeighbor;
    
    boundaryTerm->psi = &psi;
    
    nv1 = msh.p[i].size();
    // loop over all edges of the polygon we're in
    for (j = 0; j < nv1; j++)
    {
        a1 = msh.p[i][j];
        b1 = msh.p[i][(j+1)%nv1];
        
        neighbor = i;
        
        // loop over all neighboring polygons to find the one that shares this edge
        for (k = 0; k < msh.p2p[i].size(); k++)
        {
            i2 = msh.p2p[i][k];
            nv2 = msh.p[i2].size();
            // loop over all edges of the neighboring polygon
            for (l = 0; l < nv2; l++)
            {
                a2 = msh.p[i2][l];
                b2 = msh.p[i2][(l+1)%nv2];
                // have we found a match?
                if ((a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2))
                {
                    neighbor = i2;
                    break;
                }
            }
            
            // if we've found a match, don't need to keep looking
            if (neighbor != i) break;
        }
        
        msh.getOutwardNormal(i, a1, b1, boundaryTerm->nx, boundaryTerm->ny);
        boundaryTerm->iPlus = neighbor;
        boundaryTerm->UPlus = U.a.slice(neighbor);
        
        integ += msh.lineIntegral(*boundaryTerm, a1, b1);
    }
    
    return integ;
}

vec Equation::volumeIntegral(int i, const vec &psi_x, const vec &psi_y)
{
    volumeTerm->psi_x = &psi_x;
    volumeTerm->psi_y = &psi_y;
    
    return msh.polygonIntegral(*volumeTerm, i);
}

MeshFn Equation::assemble(const MeshFn &u, double t/* = 0 */)
{
    int i, j;
    int deg = u.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    double w, h;
    vec psi = zeros<vec>(basisSize);
    vec psi_x;
    vec psi_y;
    
    MeshFn b(msh, deg, nc);
    
    volumeTerm->m = deg+1;
    boundaryTerm->m = deg+1;
    
    // loop over all polygons
    for (i = 0; i < msh.np; i++)
    {
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        
        volumeTerm->i = i;
        volumeTerm->U = u.a.slice(i);
        boundaryTerm->UMinus = u.a.slice(i);
        boundaryTerm->iMinus = i;
        
        for (j = 0; j < basisSize; j++)
        {
            psi[j] = 1.0;
            
            psi_x = LegDerX(deg + 1, psi) * 2.0/w;
            psi_y = LegDerY(deg + 1, psi) * 2.0/h;
            
            b.a.slice(i).row(j) = volumeIntegral(i, psi_x, psi_y).t();
            b.a.slice(i).row(j) -= boundaryIntegral(i, psi, u).t();
          
            psi[j] = 0.0;
        }
    }
    
    return b;
}

Jacobian Equation::jacobian(const MeshFn &f, double t)
{
    int i, j, k, l, e, i2;
    int diagonalBlock, blockIdx, neighbor;
    int deg = f.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    int nv1, nv2, a1, b1, a2, b2;
    double w, h;
    vec phi = zeros<vec>(basisSize);
    vec psi = zeros<vec>(basisSize);
    vec psi_x;
    vec psi_y;
    uvec componentIndices = linspace<uvec>(0, (nc-1)*basisSize, nc);
    
    volumeJacobian->m = deg+1;
    boundaryDerivative->m = deg+1;
    
    Jacobian J(msh, f.deg, nc);
    
    for (i = 0; i < msh.np; i++)
    {
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        diagonalBlock = J.rowBlock[i];
        
        volumeJacobian->i = i;
        volumeJacobian->U = f.a.slice(i);
        boundaryDerivative->iPsi = i;
        boundaryDerivative->U = f.a.slice(i);
        
        for (j = 0; j < basisSize; j++)
        {
            phi(j) = 1;
            
            volumeJacobian->phi = &phi;
            boundaryDerivative->phi = &phi;
            
            for (k = 0; k < basisSize; k++)
            {
                psi(k) = 1;
                psi_x = LegDerX(deg + 1, psi) * 2.0/w;
                psi_y = LegDerY(deg + 1, psi) * 2.0/h;
                
                volumeJacobian->psi_x = &psi_x;
                volumeJacobian->psi_y = &psi_y;
                boundaryDerivative->psi = &psi;
                
                J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) 
                    += msh.polygonIntegral(*volumeJacobian, i);
                
                nv1 = msh.p[i].size();
    
                // loop over all edges of the polygon we're in
                for (e = 0; e < nv1; e++)
                {
                    a1 = msh.p[i][e];
                    b1 = msh.p[i][(e+1)%nv1];
                    
                    // get the outward facing normal vector
                    msh.getOutwardNormal(i, a1, b1, boundaryDerivative->nx,
                                         boundaryDerivative->ny);
                    
                    // find the neighboring polygon (or otherwise we're on the exterior)
                    neighbor = i;
                    boundaryDerivative->neighbor = i;
                    for (blockIdx = J.rowBlock[i] + 1; blockIdx < J.rowBlock[i+1]; blockIdx++)
                    {
                        i2 = J.colIndices[blockIdx];
                        
                        nv2 = msh.p[i2].size();
                        
                        // loop over all edges of the neighboring polygon
                        for (l = 0; l < nv2; l++)
                        {
                            a2 = msh.p[i2][l];
                            b2 = msh.p[i2][(l+1)%nv2];
                            // have we found a match?
                            if ((a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2))
                            {
                                neighbor = i2;
                                break;
                            }
                        }
                        
                        // compute the contribution of the neighbor on this element
                        if (neighbor != i)
                        {
                            boundaryDerivative->neighbor = neighbor;
                            boundaryDerivative->UNeighbor = f.a.slice(neighbor);
                            boundaryDerivative->iPhi = neighbor;
                            J.blocks[blockIdx](componentIndices+k, componentIndices+j) -= 
                                msh.lineIntegral(*boundaryDerivative, a1, b1);
                            
                            break;
                        }
                    }
                    
                    boundaryDerivative->iPhi = i;
                    J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) -= 
                        msh.lineIntegral(*boundaryDerivative, a1, b1);
                    
                }
                
                psi(k) = 0;
            }
            phi(j) = 0;
        }
    }
    
    return J;
}
