#include <iostream>
#include "Advection.h"
#include "Legendre.h"

using namespace arma;

double beta_x(double x, double y)
{
    return 2*y - 1;
}

double beta_y(double x, double y)
{
    return -2*x + 1;
}

double Advection::betaUDotGradPsi::operator()(double x, double y) const
{
    double xx, yy;
    msh.getLocalCoordinates(i, x, y, xx, yy);
    
    return Leg2D(xx, yy, m, u)*(beta_x(x, y)*Leg2D(xx, yy, m, *psi_x) 
                              + beta_y(x, y)*Leg2D(xx, yy, m, *psi_y));
}

double Advection::uPsiBetaDotN::operator()(double x, double y) const
{
    double xMinus, xPlus, yMinus, yPlus;
    double betaDotN;
    double psiVal;
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    
    psiVal = Leg2D(xMinus, yMinus, m, *psi);
    betaDotN = (nx*beta_x(x, y) + ny*beta_y(x, y));
    
    if (betaDotN > 0)
    {
        return betaDotN*psiVal*Leg2D(xMinus, yMinus, m, uMinus);
    } else
    {
        if (iMinus == iPlus) return 0;
        
        msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
        return betaDotN*psiVal*Leg2D(xPlus, yPlus, m, uPlus);
    }
}

Advection::Advection(PolyMesh &m)
: Equation(m), volumeTerm(betaUDotGradPsi(m)), boundaryTerm(uPsiBetaDotN(m))
{ }

double Advection::volumeIntegral(int i, vec &psi_x, vec &psi_y)
{
    volumeTerm.psi_x = &psi_x;
    volumeTerm.psi_y = &psi_y;
    
    return msh.polygonIntegral(volumeTerm, i);
}

double Advection::boundaryIntegral(int i, vec &psi, const MeshFn &u)
{
    int a1, b1, a2, b2;
    int nv1, nv2;
    int j, k, l, i2;
    int neighbor;
    double integ;
    
    integ = 0;
    
    boundaryTerm.psi = &psi;
    
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
        
        msh.getOutwardNormal(i, a1, b1, boundaryTerm.nx, boundaryTerm.ny);
        boundaryTerm.iPlus = neighbor;
        boundaryTerm.uPlus = u.a.tube(neighbor, 0);
        
        integ += msh.lineIntegral(boundaryTerm, a1, b1);
    }
    
    return integ;
}

MeshFn Advection::assemble(const MeshFn &u, double t/* = 0 */)
{
    int i, j;
    int deg = u.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    double w, h;
    vec psi = zeros<vec>(basisSize);
    vec psi_x;
    vec psi_y;
    
    // only one component for advection equation
    MeshFn b(msh, deg, 1);
    
    volumeTerm.m = deg+1;
    boundaryTerm.m = deg+1;
    
    // loop over all polygons
    for (i = 0; i < msh.np; i++)
    {
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        
        volumeTerm.i = i;
        volumeTerm.u = u.a.tube(i, 0);
        boundaryTerm.uMinus = u.a.tube(i, 0);
        boundaryTerm.iMinus = i;
        
        for (j = 0; j < basisSize; j++)
        {
            psi[j] = 1.0;
            
            psi_x = LegDerX(deg + 1, psi) * 2.0/w;
            psi_y = LegDerY(deg + 1, psi) * 2.0/h;
            
            b.a(i, 0, j) = volumeIntegral(i, psi_x, psi_y);
            b.a(i, 0, j) -= boundaryIntegral(i, psi, u);
          
            psi[j] = 0.0;
        }
    }
    
    return b;
}
