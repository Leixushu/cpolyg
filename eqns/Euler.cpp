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

mat Euler::FluxDotGradPsi::operator()(double x, double y) const
{
    double xx, yy, c, u, v;
    EulerVariables vars, flux_x, flux_y;
    
    msh.getLocalCoordinates(i, x, y, xx, yy);
    
    vars = computeVariables(xx, yy, m, U);
    
    flux(vars, gamma, flux_x, flux_y, c, u, v);
    
    return Leg2D(xx, yy, m, *psi_x)*flux_x
         + Leg2D(xx, yy, m, *psi_y)*flux_y;
}

mat Euler::LaxFriedrichsFlux::operator()(double x, double y) const
{
    EulerVariables varsMinus, varsPlus, flux_xMinus, flux_yMinus, flux_xPlus, flux_yPlus;
    double xMinus, yMinus, xPlus, yPlus;
    double cMinus, uMinus, vMinus, cPlus, uPlus, vPlus;
    double rhoPlus, rhoEPlus;
    double psiVal;
    double VDotNMinus, VDotNPlus;
    double alpha;
    
    if (iPlus != iMinus || exact == NULL)
    {
        msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
        varsPlus = computeVariables(xPlus, yPlus, m, UPlus);
    } else
    {
        rhoPlus = exact->rho(x, y);
        uPlus = exact->u(x, y);
        vPlus = exact->v(x, y);
        rhoEPlus = exact->rhoE(x, y);
        
        varsPlus[0] = rhoPlus;
        varsPlus[1] = rhoPlus*uPlus;
        varsPlus[2] = rhoPlus*vPlus;
        varsPlus[3] = rhoEPlus;
    }
    
    flux(varsPlus, gamma, flux_xPlus, flux_yPlus, cPlus, uPlus, vPlus);
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    psiVal = Leg2D(xMinus, yMinus, m, *psi);
    varsMinus = computeVariables(xMinus, yMinus, m, UMinus);
    flux(varsMinus, gamma, flux_xMinus, flux_yMinus, cMinus, uMinus, vMinus);
    
    VDotNMinus = uMinus*nx + vMinus*ny;
    VDotNPlus = uPlus*nx + vPlus*ny;
    
    alpha = max({fabs(VDotNMinus - cMinus), fabs(VDotNMinus), fabs(VDotNMinus + cMinus),
                 fabs(VDotNPlus - cPlus), fabs(VDotNPlus), fabs(VDotNPlus + cPlus)});
    
    return 0.5*((flux_xPlus + flux_xMinus)*nx + (flux_yPlus + flux_yMinus)*ny
                - alpha*(varsPlus - varsMinus))*psiVal;
}

Euler::Euler(PolyMesh &a_msh, double a_gamma) : Equation(a_msh)
{
    nc = kEulerComponents;
    gamma = a_gamma;
    exact = NULL;
    volumeTerm = new FluxDotGradPsi(a_msh, gamma);
    boundaryTerm = new LaxFriedrichsFlux(a_msh, gamma);
}

Euler::~Euler()
{
    delete volumeTerm;
    delete boundaryTerm;
    if(exact) delete exact;
}
