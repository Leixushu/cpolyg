#include <iostream>
#include "Advection.h"
#include "Legendre.h"
#include "Timer/CH_Timer.H"

using namespace arma;

double Advection::betaX(double x, double y) const
{
    return 2*y - 1;
}

double Advection::betaY(double x, double y) const
{
    return -2*x + 1;
}

mat Advection::computeVolumeTerm(double x, double y)
{
    double xx, yy;
    vec::fixed<1> result;
    msh.getLocalCoordinates(iMinus, x, y, xx, yy);
    
    result(0) = Leg2D(xx, yy, m, UMinus)*(betaX(x, y)*Leg2D(xx, yy, m, psi_x) 
                                        + betaY(x, y)*Leg2D(xx, yy, m, psi_y));
    return result;
}

mat Advection::computeBoundaryTerm(double x, double y)
{
    double xMinus, xPlus, yMinus, yPlus;
    double betaDotN;
    double psiVal;
    vec::fixed<1> result;
    
    msh.getLocalCoordinates(iMinus, x, y, xMinus, yMinus);
    
    psiVal = Leg2D(xMinus, yMinus, m, psi);
    betaDotN = (nx*betaX(x, y) + ny*betaY(x, y));
    
    if (betaDotN > 0)
    {
        result(0) = betaDotN*psiVal*Leg2D(xMinus, yMinus, m, UMinus);
        return result;
    } else
    {
        // negative iPlus indicates exterior boundary
        if (iPlus < 0)
        {
            result = bc.boundaryValue(x, y, UPlus, iPlus);
            result(0) *= betaDotN*psiVal;
        } else
        {
            msh.getLocalCoordinates(iPlus, x, y, xPlus, yPlus);
            result(0) = betaDotN*psiVal*Leg2D(xPlus, yPlus, m, UPlus);
        }
        return result;
    }
}

mat Advection::computeVolumeJacobian(double x, double y)
{
    double xx, yy;
    vec::fixed<1> result;
    msh.getLocalCoordinates(iPsi, x, y, xx, yy);
    
    result(0) = Leg2D(xx, yy, m, phi)*(betaX(x, y)*Leg2D(xx, yy, m, psi_x) 
                                     + betaY(x, y)*Leg2D(xx, yy, m, psi_y));
    return result;
}

mat Advection::computeBoundaryJacobian(double x, double y)
{
    double xPhi, yPhi, xPsi, yPsi;
    double betaDotN;
    vec::fixed<1> result;
    int sgn;
    
    if (iPhi == iPsi) sgn = 1; else sgn = -1;
    
    betaDotN = nx*betaX(x, y) + ny*betaY(x, y);
    
    if (sgn*betaDotN > 0)
    {
        msh.getLocalCoordinates(iPsi, x, y, xPsi, yPsi);
        if (iPhi >= 0)
        {
            msh.getLocalCoordinates(iPhi, x, y, xPhi, yPhi);
            result(0) = betaDotN*Leg2D(xPsi, yPsi, m, psi)*Leg2D(xPhi, yPhi, m, phi);
        } else
        {
            result = bc.boundaryValue(x, y, phi, iPhi);
            result(0) *= betaDotN*Leg2D(xPsi, yPsi, m, psi);
        }
    } else
    {
        result(0) = 0;
    }
    
    return result;
}
