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

void Euler::flux(const EulerVariables &vars, vec &flux_x, vec &flux_y, 
                 double &c, double &u, double &v)
{
    double rho, rhoU, rhoV, rhoE, P;
    
    rho = vars(0);
    rhoU = vars(1);
    rhoV = vars(2);
    rhoE = vars(3);
    
    u = rhoU/rho;
    v = rhoV/rho;
    
    P = (gamma - 1)*(rhoE - 0.5*(rhoU*rhoU + rhoV*rhoV)/rho);
    
    flux_x[0] = rhoU;
    flux_x[1] = rhoU*u + P;
    flux_x[2] = rhoU*v;
    flux_x[3] = (rhoE + P)*u;
    
    flux_y[0] = rhoV;
    flux_y[1] = rhoU*v;
    flux_y[2] = rhoV*v + P;
    flux_y[3] = (rhoE + P)*v; 
    
    c = sqrt(gamma*P/rho);
}

void Euler::fluxJacobian(const EulerVariables &vars, EulerJacobian &J1, EulerJacobian &J2)
{
    double rho, rhoU, rhoV, rhoE, u, v, E, P, usq, e, H;
    
    rho = vars(0);
    rhoU = vars(1);
    rhoV = vars(2);
    rhoE = vars(3);
    
    u = rhoU/rho;
    v = rhoV/rho;
    E = rhoE/rho;
    
    usq = u*u + v*v;
    e = E - 0.5*usq;
    P = (gamma - 1)*rho*e;
    H = E + P/rho;
    
    J1 << 0                         << 1                 << 0              << 0       << endr
       << -u*u + 0.5*(gamma-1)*usq  << (3-gamma)*u       << -(gamma-1)*v   << gamma-1 << endr
       << -u*v                      << v                 << u              << 0       << endr
       << u*(0.5*(gamma-1)*usq - H) << H - (gamma-1)*usq << -(gamma-1)*u*v << gamma*u << endr;
       
    
    J2 << 0                           << 0                << 1                   << 0         << endr
       << -u*v                        << v                << u                   << 0         << endr
       << -v*v + 0.5*(gamma - 1)*usq  << -(gamma - 1)*u   << (3 - gamma)*v       << gamma - 1 << endr
       << v*(0.5*(gamma - 1)*usq - H) << -(gamma - 1)*u*v << H - (gamma - 1)*v*v << gamma*v   << endr;
}

EulerVariables Euler::alphaDerivative(const EulerVariables &vars1, const EulerVariables &vars2, 
                                      double &alpha)
{
    double u1, v1, u2, v2, c1, c2, P1, P2, vDotN1, vDotN2, rho, usq, a, b, E;
    uword alphaIdx;
    vec::fixed<6> eigs;
    EulerVariables derivative;
    derivative.zeros();
    
    u1 = vars1(1)/vars1(0);
    v1 = vars1(2)/vars1(0);
    u2 = vars2(1)/vars2(0);
    v2 = vars2(2)/vars2(0);
    
    P1 = (gamma - 1)*(vars1(3) - 0.5*(vars1(1)*vars1(1) + vars1(2)*vars1(2))/vars1(0));
    c1 = sqrt(gamma*P1/vars1(0));
    
    P2 = (gamma - 1)*(vars2(3) - 0.5*(vars2(1)*vars2(1) + vars2(2)*vars2(2))/vars2(0));
    c2 = sqrt(gamma*P2/vars2(0));
    
    vDotN1 = u1*nx + v1*ny;
    vDotN2 = u2*nx + v2*ny;
    
    eigs = {vDotN1 - c1, vDotN1, vDotN1 + c1,
            vDotN2 - c2, vDotN2, vDotN2 + c2};
    
    alpha = abs(eigs).max(alphaIdx);
    
    if (alphaIdx <= 2)
    {
        rho = vars1(0);
        E = vars1(3)/rho;
        
        derivative(0) = -u1/rho*nx - v1/rho*ny;
        derivative(1) = 1/rho*nx;
        derivative(2) = 1/rho*ny;
        
        usq = u1*u1 + v1*v1;
        
        a = sqrt(gamma*gamma - gamma);
        b = rho*sqrt(4*E - 2*usq);
        
        if (alphaIdx == 0)
        {
            derivative(0) -= a*(usq - E)/b;
            derivative(1) -= -a*u1/b;
            derivative(2) -= -a*v1/b;
            derivative(3) -= a/b;
        } else if (alphaIdx == 2)
        {
            derivative(0) += a*(usq - E)/b;
            derivative(1) += -a*u1/b;
            derivative(2) += -a*v1/b;
            derivative(3) += a/b;
        }
        
        if (eigs(alphaIdx) < 0)
        {
            derivative *= -1;
        }
    }
    
    return derivative;
}

mat Euler::computeVolumeTerm(double x, double y)
{
    double xx, yy, c, u, v;
    EulerVariables vars, flux_x, flux_y;
    
    msh.getLocalCoordinates(iMinus, x, y, xx, yy);
    
    vars = computeVariables(xx, yy, m, UMinus);
    
    flux(vars, flux_x, flux_y, c, u, v);
    
    return Leg2D(xx, yy, m, psi_x)*flux_x
         + Leg2D(xx, yy, m, psi_y)*flux_y;
}

mat Euler::computeBoundaryTerm(double x, double y)
{
    EulerVariables varsMinus, varsPlus, flux_xMinus, flux_yMinus, flux_xPlus, flux_yPlus;
    double xMinus, yMinus, xPlus, yPlus;
    double cMinus, uMinus, vMinus, cPlus, uPlus, vPlus;
    double psiVal;
    double VDotNMinus, VDotNPlus;
    double alpha;
    
    if (iPlus >= 0)
    {
        msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
        varsPlus = computeVariables(xPlus, yPlus, m, UPlus);
    } else
    {
        varsPlus = bc.boundaryValue(x, y, UPlus, iPlus);
    }
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    psiVal = Leg2D(xMinus, yMinus, m, psi);
    varsMinus = computeVariables(xMinus, yMinus, m, UMinus);
    
    flux(varsMinus, flux_xMinus, flux_yMinus, cMinus, uMinus, vMinus);
    flux(varsPlus, flux_xPlus, flux_yPlus, cPlus, uPlus, vPlus);
    
    VDotNMinus = uMinus*nx + vMinus*ny;
    VDotNPlus = uPlus*nx + vPlus*ny;
    
    alpha = max({fabs(VDotNMinus - cMinus), fabs(VDotNMinus), fabs(VDotNMinus + cMinus),
                 fabs(VDotNPlus - cPlus), fabs(VDotNPlus), fabs(VDotNPlus + cPlus)});
    
    return 0.5*((flux_xPlus + flux_xMinus)*nx + (flux_yPlus + flux_yMinus)*ny
                + alpha*(varsMinus - varsPlus))*psiVal;
}

mat Euler::computeVolumeJacobian(double x, double y)
{
    double xx, yy, phiVal;
    EulerVariables vars;
    EulerJacobian J1, J2;
    
    msh.getLocalCoordinates(iPsi, x, y, xx, yy);
    
    vars = computeVariables(xx, yy, m, U);
    
    fluxJacobian(vars, J1, J2);
    
    phiVal = Leg2D(xx, yy, m, phi);
    return phiVal*(Leg2D(xx, yy, m, psi_x)*J1 + Leg2D(xx, yy, m, psi_y)*J2);
}

mat Euler::computeBoundaryJacobian(double x, double y)
{
    double xPhi, yPhi, xPsi, yPsi, xNeighbor, yNeighbor, psiVal, phiVal;
    EulerVariables vars, vars2, alphaPrime;
    EulerJacobian J1, J2;
    int sgn;
    double alpha;
    
    msh.getLocalCoordinates(iPsi, x, y, xPsi, yPsi);
    psiVal = Leg2D(xPsi, yPsi, m, psi);
    vars = computeVariables(xPsi, yPsi, m, U);
    
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
        vars2 = computeVariables(xNeighbor, yNeighbor, m, UNeighbor);
    } else
    {
        vars2 = bc.boundaryValue(x, y, UNeighbor, neighbor);
    }
    
    if (iPsi == iPhi)
    {
        fluxJacobian(vars, J1, J2);
        sgn = 1;
    } else
    {
        fluxJacobian(vars2, J1, J2);
        sgn = -1;
    }
    
    alphaPrime = alphaDerivative(vars, vars2, alpha);
    
    return 0.5*phiVal*psiVal*(J1*nx + J2*ny + sgn*alpha*Id
                            + (vars - vars2)*alphaPrime.t());
}

Euler::Euler(PolyMesh &a_msh, BoundaryConditions a_bc, double g)
: Equation(a_msh, a_bc, kEulerComponents), gamma(g)
{ }
