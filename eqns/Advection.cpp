#include <iostream>
#include "Advection.h"
#include "Legendre.h"
#include "Timer/CH_Timer.H"

using namespace arma;

double beta_x(double x, double y)
{
    //return 2;
    //return 2*y - 1;
    return 0.866025;
}

double beta_y(double x, double y)
{
    //return 1;
    //return -2*x + 1;
    return 0.5;
}

mat Advection::betaUDotGradPsi::operator()(double x, double y) const
{
    double xx, yy;
    vec::fixed<1> result;
    msh.getLocalCoordinates(i, x, y, xx, yy);
    
    result(0) = Leg2D(xx, yy, m, U)*(beta_x(x, y)*Leg2D(xx, yy, m, *psi_x) 
                                   + beta_y(x, y)*Leg2D(xx, yy, m, *psi_y));
    return result;
}

mat Advection::JacobianBetaUDotGradPsi::operator()(double x, double y) const
{
    double xx, yy;
    vec::fixed<1> result;
    msh.getLocalCoordinates(i, x, y, xx, yy);
    
    result(0) = Leg2D(xx, yy, m, *phi)*(beta_x(x, y)*Leg2D(xx, yy, m, *psi_x) 
                                      + beta_y(x, y)*Leg2D(xx, yy, m, *psi_y));
    return result;
}

mat Advection::uPsiBetaDotN::operator()(double x, double y) const
{
    double xMinus, xPlus, yMinus, yPlus;
    double betaDotN;
    double psiVal;
    vec::fixed<1> result;
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    
    psiVal = Leg2D(xMinus, yMinus, m, *psi);
    betaDotN = (nx*beta_x(x, y) + ny*beta_y(x, y));
    
    if (betaDotN > 0)
    {
        result(0) = betaDotN*psiVal*Leg2D(xMinus, yMinus, m, UMinus);
        return result;
    } else
    {
        if (iMinus == iPlus)
        {
            result(0) = 0;
        } else
        {
            msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
            result(0) = betaDotN*psiVal*Leg2D(xPlus, yPlus, m, UPlus);
        }
        return result;
    }
}

mat Advection::phiPsiBetaDotN::operator()(double x, double y) const
{
    double xPhi, yPhi, xPsi, yPsi;
    double betaDotN;
    vec::fixed<1> result;
    int sgn;
    
    if (iPhi == iPsi) sgn = 1; else sgn = -1;
    
    betaDotN = nx*beta_x(x, y) + ny*beta_y(x, y);
    
    if (sgn*betaDotN > 0)
    {
        msh.getLocalCoordinates(iPhi, x, y, xPhi, yPhi);
        msh.getLocalCoordinates(iPsi, x, y, xPsi, yPsi);
        
        result(0) = betaDotN*Leg2D(xPsi, yPsi, m, *psi)*Leg2D(xPhi, yPhi, m, *phi);
        
    } else
    {
        result(0) = 0;
    }
    
    return result;
}

Advection::Advection(PolyMesh &a_msh)
: Equation(a_msh)
{
    nc = 1;
    volumeTerm = new betaUDotGradPsi(a_msh);
    volumeJacobian = new JacobianBetaUDotGradPsi(a_msh);
    boundaryTerm = new uPsiBetaDotN(a_msh);
    boundaryDerivative = new phiPsiBetaDotN(a_msh);
}

Advection::~Advection()
{
    delete volumeTerm;
    delete volumeJacobian;
    delete boundaryTerm;
    delete boundaryDerivative;
}
