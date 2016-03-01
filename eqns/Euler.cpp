#include "Euler.h"
#include "Legendre.h"

using namespace arma;
using namespace std;

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

mat Euler::fluxFunction(const vec &val, double x, double y)
{
    mat result(kEulerComponents, 2);
    EulerVariables flux_x, flux_y;
    double c, u, v;
    
    flux(val, flux_x, flux_y, c, u, v);
    result.col(0) = flux_x;
    result.col(1) = flux_y;
    
    return result;
}

vec Euler::numericalFluxFunction(const vec &varsMinus, const vec &varsPlus, 
    double x, double y, double nx, double ny)
{
    EulerVariables flux_xMinus, flux_yMinus, flux_xPlus, flux_yPlus;
    double cMinus, uMinus, vMinus, cPlus, uPlus, vPlus;
    double VDotNMinus, VDotNPlus;
    double alpha;
    
    flux(varsMinus, flux_xMinus, flux_yMinus, cMinus, uMinus, vMinus);
    flux(varsPlus, flux_xPlus, flux_yPlus, cPlus, uPlus, vPlus);
    
    VDotNMinus = uMinus*nx + vMinus*ny;
    VDotNPlus = uPlus*nx + vPlus*ny;
    
    alpha = max({fabs(VDotNMinus - cMinus), fabs(VDotNMinus), fabs(VDotNMinus + cMinus),
                 fabs(VDotNPlus - cPlus), fabs(VDotNPlus), fabs(VDotNPlus + cPlus)});
    
    return 0.5*((flux_xPlus + flux_xMinus)*nx + (flux_yPlus + flux_yMinus)*ny
                + alpha*(varsMinus - varsPlus));
}

cube Euler::fluxJacobian(const vec &vars, double x, double y)
{
    EulerJacobian J1, J2;
    cube::fixed<kEulerComponents, kEulerComponents, 2> result;
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
    
    result.slice(0) =
        {{0, 1, 0, 0},
         {-u*u + 0.5*(gamma-1)*usq, (3-gamma)*u, -(gamma-1)*v, gamma-1},
         {-u*v, v, u, 0},
         {u*(0.5*(gamma-1)*usq - H), H - (gamma-1)*usq, -(gamma-1)*u*v, gamma*u}};
    result.slice(1) = 
        {{0, 0, 1, 0},
         {-u*v, v, u, 0},
         {-v*v + 0.5*(gamma - 1)*usq, -(gamma - 1)*u, (3 - gamma)*v, gamma - 1},
         {v*(0.5*(gamma - 1)*usq - H), -(gamma - 1)*u*v, H - (gamma - 1)*v*v, gamma*v}};
    
    return result;
}

mat Euler::numericalFluxJacobian(const vec &vars, const vec &vars2, double x, 
    double y, double nx, double ny, int sgn)
{
    EulerVariables alphaPrime;
    double alpha;
    
    cube J = fluxJacobian(vars, x, y);
    alphaPrime = alphaDerivative(vars, vars2, alpha);
    
    return 0.5*(J.slice(0)*nx + J.slice(1)*ny + sgn*alpha*Id
                + (vars - vars2)*alphaPrime.t());
}
