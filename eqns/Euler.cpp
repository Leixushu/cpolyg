#include "Euler.h"
#include "Legendre.h"

using namespace arma;
using namespace std;

EulerVariables Euler::computeVariables(double x, double y, int m, const mat &U)
{
    EulerVariables values;
    values(0) = Leg2D(x, y, m, U.col(0));
    values(1) = Leg2D(x, y, m, U.col(1));
    values(2) = Leg2D(x, y, m, U.col(2));
    values(3) = Leg2D(x, y, m, U.col(3));
    
    return values;
}

void Euler::flux(const EulerVariables &U, double gamma, vec &flux_x, vec &flux_y, 
                 double &c, double &u, double &v)
{
    double rho, rhoU, rhoV, rhoE, P;
    
    rho = U(0);
    rhoU = U(1);
    rhoV = U(2);
    rhoE = U(3);
    
    u = rhoU/rho;
    v = rhoV/rho;
    
    P = (gamma - 1)*(rhoE - 0.5*(rhoU*rhoU + rhoV*rhoV)/rho);
    
    flux_x[0] = rhoU;
    flux_x[1] = rhoU*u + P;
    flux_x[2] = rhoU*v;
    flux_x[3] = rhoE*u + P*u;
    
    flux_y[0] = rhoV;
    flux_y[1] = rhoU*v;
    flux_y[2] = rhoV*v + P;
    flux_y[3] = rhoE*v + P*v; 
    
    c = sqrt(gamma*P/rho);
}

vec Euler::FluxDotGradPsi::operator()(double x, double y) const
{
    double xx, yy, c, u, v;
    EulerVariables vars, flux_x, flux_y;
    
    msh.getLocalCoordinates(i, x, y, xx, yy);
    
    vars = computeVariables(xx, yy, m, *U);
    
    flux(vars, gamma, flux_x, flux_y, c, u, v);
    
    return Leg2D(xx, yy, m, *psi_x)*flux_x
         + Leg2D(xx, yy, m, *psi_y)*flux_y;
}

vec Euler::LaxFriedrichsFlux::operator()(double x, double y) const
{
    EulerVariables varsMinus, varsPlus, flux_xMinus, flux_yMinus, flux_xPlus, flux_yPlus;
    double xMinus, yMinus, xPlus, yPlus;
    double cMinus, uMinus, vMinus, cPlus, uPlus, vPlus;
    double psiVal;
    double VDotNMinus, VDotNPlus;
    double alpha;
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
    
    psiVal = Leg2D(xMinus, yMinus, m, *psi);
    
    varsMinus = computeVariables(xMinus, yMinus, m, *UMinus);
    varsPlus = computeVariables(xPlus, yPlus, m, *UPlus);
    
    flux(varsMinus, gamma, flux_xMinus, flux_yMinus, cMinus, uMinus, vMinus);
    flux(varsPlus, gamma, flux_xPlus, flux_yPlus, cPlus, uPlus, vPlus);
    
    VDotNMinus = uMinus*nx + vMinus*ny;
    VDotNPlus = uPlus*nx + vPlus*ny;
    
    alpha = max({fabs(VDotNMinus - cMinus), fabs(VDotNMinus), fabs(VDotNMinus + cMinus),
                 fabs(VDotNPlus - cPlus), fabs(VDotNPlus), fabs(VDotNPlus + cPlus)});
    
    return 0.5*((flux_xPlus + flux_xMinus)*nx + (flux_yPlus + flux_yMinus)*ny
                - alpha*(varsPlus - varsMinus));
}

vec Euler::boundaryIntegral(int i, const vec &psi, const MeshFn &U)
{
    int a1, b1, a2, b2;
    int nv1, nv2;
    int j, k, l, i2;
    int neighbor;
    vec::fixed<kEulerComponents> integ;
    mat UNeighbor;
    
    integ.fill(0);
    
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
        
        UNeighbor = U.a.tube(neighbor, 0, neighbor, kEulerComponents-1);
        
        msh.getOutwardNormal(i, a1, b1, boundaryTerm.nx, boundaryTerm.ny);
        boundaryTerm.iPlus = neighbor;
        boundaryTerm.UPlus = &UNeighbor;
        
        integ += msh.lineIntegral(boundaryTerm, a1, b1);
    }
    
    return integ;
}

vec Euler::volumeIntegral(int i, const vec &psi_x, const vec &psi_y)
{
    volumeTerm.psi_x = &psi_x;
    volumeTerm.psi_y = &psi_y;
    
    return msh.polygonIntegral(volumeTerm, i);
}

Euler::Euler(PolyMesh &m, double g)
: Equation(m), gamma(g), volumeTerm(m, g), boundaryTerm(m, g)
{ };

MeshFn Euler::assemble(const MeshFn &U)
{
    int i, j;
    int deg = U.deg;
    int basisSize = (deg+1)*(deg+2)/2;
    double w, h;
    mat ULocal;
    
    vec psi = zeros<vec>(basisSize);
    vec psi_x;
    vec psi_y;
    
    // four components for 2D Euler equations
    MeshFn b(msh, deg, kEulerComponents);
    
    volumeTerm.m = deg+1;
    boundaryTerm.m = deg+1;
    
    // loop over all polygons
    for (i = 0; i < msh.np; i++)
    {
        w = msh.bb[i][2];
        h = msh.bb[i][3];
        
        // extract the local coefficients for this component
        ULocal = U.a.tube(i, 0, i, kEulerComponents-1);
        
        volumeTerm.i = i;
        volumeTerm.U = &ULocal;
        boundaryTerm.UMinus = &ULocal;
        boundaryTerm.iMinus = i;
        
        for (j = 0; j < basisSize; j++)
        {
            psi[j] = 1.0;
            
            psi_x = LegDerX(deg + 1, psi) * 2.0/w;
            psi_y = LegDerY(deg + 1, psi) * 2.0/h;
            
            b.a.subcube(i, 0, j, i, kEulerComponents-1, j) 
                = volumeIntegral(i, psi, ULocal);
            b.a.subcube(i, 0, j, i, kEulerComponents-1, j) 
                -= boundaryIntegral(i, psi, U);
          
            psi[j] = 0.0;
        }
    }
    
    return b;

}
