#include "Equation.h"
#include "Quadrature.h"

using namespace arma;
using namespace std;

Equation::Equation(PolyMesh &a_msh) : msh(a_msh)
{
    volumeTerm = NULL;
    boundaryTerm = NULL;
    volumeJacobian = NULL;
    boundaryDerivative = NULL;
}

vec Equation::boundaryIntegral(int i, const vec &psi, const MeshFn &U, int deg)
{
    int nv1, a1, b1;
    int j;
    vec integ = zeros<vec>(nc);
    
    boundaryTerm->psi = &psi;
    
    nv1 = msh.p[i].size();
    // loop over all edges of the polygon we're in
    for (j = 0; j < nv1; j++)
    {
        a1 = msh.p[i][j];
        b1 = msh.p[i][(j+1)%nv1];
        
        
        boundaryTerm->iPlus = msh.p2p[i][j];
        if(msh.p2p[i][j] >= 0)
        {
            boundaryTerm->UPlus = U.a.slice(msh.p2p[i][j]);
        } else
        {
            // how to deal with BCs in a good way???
            //int p2 = dynamic_cast<PeriodicMesh &>(msh).bc[-msh.p2p[i][j]-1].p2;
            //boundaryTerm->UPlus = U.a.slice(p2);
        }
        
        msh.getOutwardNormal(i, a1, b1, boundaryTerm->nx, boundaryTerm->ny);
        integ += Quadrature::lineIntegral(msh, *boundaryTerm, a1, b1, deg*2);
    }
    
    return integ;
}

vec Equation::volumeIntegral(int i, const vec &psi_x, const vec &psi_y, int deg)
{
    volumeTerm->psi_x = &psi_x;
    volumeTerm->psi_y = &psi_y;
    
    return Quadrature::polygonIntegral(msh, *volumeTerm, i, deg*2);
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
            
            b.a.slice(i).row(j) = volumeIntegral(i, psi_x, psi_y, deg).t();
            b.a.slice(i).row(j) -= boundaryIntegral(i, psi, u, deg).t();
          
            psi[j] = 0.0;
        }
    }
    
    return b;
}

Jacobian Equation::jacobian(const MeshFn &f, double t)
{
    int i, j, k, e;
    int diagonalBlock, blockIdx;
    int deg = f.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    int nv1, a1, b1, neighbor;
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
                    += Quadrature::polygonIntegral(msh, *volumeJacobian, i, deg*2);
                
                nv1 = msh.p[i].size();
                
                blockIdx = diagonalBlock + 1;
                
                // loop over all edges of the polygon we're in
                for (e = 0; e < nv1; e++)
                {
                    a1 = msh.p[i][e];
                    b1 = msh.p[i][(e+1)%nv1];
                    
                    // get the outward facing normal vector
                    msh.getOutwardNormal(i, a1, b1, boundaryDerivative->nx,
                                         boundaryDerivative->ny);
                    
                    // get the neighbor corresponding to this edge
                    neighbor = msh.p2p[i][e];
                    boundaryDerivative->neighbor = neighbor;
                    
                    // negative index indicates exterior edge
                    if (neighbor >= 0)
                    {
                        boundaryDerivative->UNeighbor = f.a.slice(neighbor);
                        boundaryDerivative->iPhi = neighbor;
                        J.blocks[blockIdx](componentIndices+k, componentIndices+j) -= 
                            Quadrature::lineIntegral(msh, *boundaryDerivative, a1, b1, deg*2);
                        blockIdx++;
                    }
                    
                    boundaryDerivative->iPhi = i;
                    J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) -= 
                        Quadrature::lineIntegral(msh, *boundaryDerivative, a1, b1, deg*2);
                    
                }
                
                psi(k) = 0;
            }
            phi(j) = 0;
        }
    }
    
    return J;
}
