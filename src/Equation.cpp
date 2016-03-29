#include "Equation.h"
#include "Quadrature.h"
#include "Timer/CH_Timer.H"

using namespace arma;
using namespace std;

vec Equation::computeVariables(const mat &coeffs, double x, double y)
{
    vec values(nc);
    
    for (size_t c = 0; c < nc; c++)
    {
        values(c) = Leg2D(x, y, m, coeffs.col(c));
    }
    
    return values;
}

mat Equation::computeVolumeTerm(double x, double y)
{
    double xx, yy;
    vec result, vars; // nc
    
    msh.getLocalCoordinates(iMinus, x, y, xx, yy);
    vars = computeVariables(UMinus, xx, yy);
    mat theFlux = fluxFunction(vars, x, y);
    
    result = (theFlux.col(0)*Leg2D(xx, yy, m, psi_x)
            + theFlux.col(1)*Leg2D(xx, yy, m, psi_y));
    
    return result;
}

mat Equation::computeBoundaryTerm(double x, double y)
{
    double xMinus, xPlus, yMinus, yPlus;
    double psiVal;
    vec varsMinus, varsPlus; // nc
    vec theFlux; // nc
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    psiVal = Leg2D(xMinus, yMinus, m, psi);
    varsMinus = computeVariables(UMinus, xMinus, yMinus);
    
    if (iPlus < 0)
    {
        varsPlus = bc.boundaryValue(x, y, UPlus, iPlus);
    } else
    {
        msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
        varsPlus = computeVariables(UPlus, xPlus, yPlus);
    }
    
    theFlux = numericalFluxFunction(varsMinus, varsPlus, x, y, nx, ny);
    return psiVal*theFlux;
}

mat Equation::computeVolumeJacobian(double x, double y)
{
    vec vars; // nc
    mat result; // (nc)x(nc)
    cube theJacobian; // (nc)x(nc)x2
    double xx, yy;
    
    msh.getLocalCoordinates(iPsi, x, y, xx, yy);
    vars = computeVariables(U, xx, yy);
    theJacobian = fluxJacobian(vars, x, y);
    
    result = Leg2D(xx, yy, m, phi)*(theJacobian.slice(0)*Leg2D(xx, yy, m, psi_x) 
                                  + theJacobian.slice(1)*Leg2D(xx, yy, m, psi_y));
    return result;
}

mat Equation::computeBoundaryJacobian(double x, double y)
{
    double xPhi, yPhi, xPsi, yPsi, xNeighbor, yNeighbor;
    double psiVal, phiVal;
    vec vars2;
    int sgn;
    
    msh.getLocalCoordinates(iPsi, x, y, xPsi, yPsi);
    psiVal = Leg2D(xPsi, yPsi, m, psi);
    mat vars = computeVariables(U, xPsi, yPsi);
    
    if (iPhi >= 0)
    {
        msh.getLocalCoordinates(iPhi, x, y, xPhi, yPhi);
        phiVal = Leg2D(xPhi, yPhi, m, phi);
    } else
    {
        // if iPhi < 0, we're at a periodic boundary
        phiVal = bc.boundaryValue(x, y, phi, iPhi)(0);
    }
    
    // if we're not on an exterior edge, take the value of the neighboring cell
    if (neighbor >= 0)
    {
        msh.getLocalCoordinates(neighbor, x, y, xNeighbor, yNeighbor);
        vars2 = computeVariables(UNeighbor, xNeighbor, yNeighbor);
    } else
    {
        vars2 = bc.boundaryValue(x, y, UNeighbor, neighbor);
    }
    
    if (iPhi == iPsi) sgn = 1;
    else sgn = -1;
    
    mat jacobian = numericalFluxJacobian(vars, vars2, x, y, nx, ny, sgn);
    
    return phiVal*psiVal*jacobian;
}

vec Equation::boundaryIntegral(int i, const MeshFn &u, int deg)
{
    int nv1, a1, b1;
    int j;
    vec integ = zeros<vec>(nc);
    
    nv1 = msh.p[i].size();
    // loop over all edges of the polygon we're in
    for (j = 0; j < nv1; j++)
    {
        a1 = msh.p[i][j];
        b1 = msh.p[i][(j+1)%nv1];
        msh.getOutwardNormal(i, a1, b1, nx, ny);
        
        // get the neighboring cell
        iPlus = msh.p2p[i][j];
        // a nonnegative number indicates a neighbor
        // a negative number indicates a boundary edge
        // (but the mesh could be periodic, have to check)
        if(iPlus >= 0)
        {
            UPlus = u.a.slice(iPlus);
        } else if (msh.bc.at(iPlus).periodic)
        {
            UPlus = u.a.slice(msh.bc.at(iPlus).p2);
        }
        
        integ += Quadrature::lineIntegral(msh, boundaryTerm, a1, b1, deg*2);
    }
    
    return integ;
}

vec Equation::volumeIntegral(int i, int deg)
{
    return Quadrature::polygonIntegral(msh, volumeTerm, i, deg*2);
}

MeshFn Equation::assemble(const MeshFn &u, double a_t/* = 0 */)
{
    CH_TIMERS("Assemble");
    int i, j;
    int deg = u.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    double w, h;
    MeshFn b(msh, deg, nc);
    
    psi = zeros<vec>(basisSize);
    // todo: get rid of these stupid 'm's
    m = deg+1;
    bc.t = a_t;
    
    // loop over all polygons
    for (i = 0; i < msh.np; i++)
    {
        iMinus = i;
        
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        
        UMinus = u.a.slice(i);
        
        for (j = 0; j < basisSize; j++)
        {
            psi[j] = 1.0;
            
            psi_x = LegDerX(deg + 1, psi) * 2.0/w;
            psi_y = LegDerY(deg + 1, psi) * 2.0/h;
            
            b.a.slice(i).row(j) = volumeIntegral(i, deg).t();
            b.a.slice(i).row(j) -= boundaryIntegral(i, u, deg).t();
          
            psi[j] = 0.0;
        }
    }
    
    return b;
}

Jacobian Equation::jacobian(const MeshFn &f, double a_t)
{
    CH_TIMERS("Jacobian");
    int i, j, k, e;
    int diagonalBlock, blockIdx;
    int deg = f.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    int nv1, a1, b1;
    double w, h;
    uvec componentIndices = linspace<uvec>(0, (nc-1)*basisSize, nc);
    
    bc.t = a_t;
    
    phi = zeros<vec>(basisSize);
    psi = zeros<vec>(basisSize);
    
    m = deg+1;
    
    Jacobian J(msh, f.deg, nc);
    
    for (i = 0; i < msh.np; i++)
    {
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        diagonalBlock = J.rowBlock[i];
        
        U = f.a.slice(i);
        iPsi = i;
        
        for (j = 0; j < basisSize; j++)
        {
            phi(j) = 1;
            
            for (k = 0; k < basisSize; k++)
            {
                psi(k) = 1;
                psi_x = LegDerX(deg + 1, psi) * 2.0/w;
                psi_y = LegDerY(deg + 1, psi) * 2.0/h;
                
                J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) 
                    += Quadrature::polygonIntegral(msh, volumeJacobian, i, deg*2);
                
                nv1 = msh.p[i].size();
                
                blockIdx = diagonalBlock + 1;
                
                // loop over all edges of the polygon we're in
                for (e = 0; e < nv1; e++)
                {
                    a1 = msh.p[i][e];
                    b1 = msh.p[i][(e+1)%nv1];
                    
                    // get the outward facing normal vector
                    msh.getOutwardNormal(i, a1, b1, nx, ny);
                    
                    // get the neighbor corresponding to this edge
                    neighbor = msh.p2p[i][e];
                    iPhi = neighbor;
                    
                    // negative index indicates exterior edge
                    if (neighbor >= 0)
                    {
                        UNeighbor = f.a.slice(neighbor);
                        J.blocks[blockIdx](componentIndices+k, componentIndices+j) -= 
                            Quadrature::lineIntegral(msh, boundaryJacobian, a1, b1, deg*2);
                        blockIdx++;
                    } else if (msh.bc.at(neighbor).periodic)
                    {
                        UNeighbor = f.a.slice(msh.bc.at(neighbor).p2);
                        J.blocks[blockIdx](componentIndices+k, componentIndices+j) -= 
                            Quadrature::lineIntegral(msh, boundaryJacobian, a1, b1, deg*2);
                        blockIdx++;
                    }
                    
                    iPhi = i;
                    J.blocks[diagonalBlock](componentIndices+k, componentIndices+j) -= 
                        Quadrature::lineIntegral(msh, boundaryJacobian, a1, b1, deg*2);
                }
                
                psi(k) = 0;
            }
            phi(j) = 0;
        }
    }
    
    return J;
}

Equation::Equation(PolyMesh &a_msh, BoundaryConditions a_bc, int a_nc)
: msh(a_msh), 
  volumeTerm(*this, &Equation::computeVolumeTerm, a_nc, 1),
  boundaryTerm(*this, &Equation::computeBoundaryTerm, a_nc, 1),
  volumeJacobian(*this, &Equation::computeVolumeJacobian, a_nc, a_nc),
  boundaryJacobian(*this, &Equation::computeBoundaryJacobian, a_nc, a_nc),
  bc(a_bc), nc(a_nc)
{
}

Equation::~Equation()
{
}
