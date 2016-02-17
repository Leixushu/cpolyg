#include <iostream>
#include "Advection.h"
#include "Legendre.h"
#include "Timer/CH_Timer.H"

using namespace arma;

double beta_x(double x, double y)
{
    return 2*y - 1;
}

double beta_y(double x, double y)
{
    return -2*x + 1;
}

mat Advection::computeVolumeTerm(double x, double y)
{
    double xx, yy;
    vec::fixed<1> result;
    msh.getLocalCoordinates(iMinus, x, y, xx, yy);
    
    result(0) = Leg2D(xx, yy, m, UMinus)*(beta_x(x, y)*Leg2D(xx, yy, m, psi_x) 
                                        + beta_y(x, y)*Leg2D(xx, yy, m, psi_y));
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
    betaDotN = (nx*beta_x(x, y) + ny*beta_y(x, y));
    
    if (betaDotN > 0)
    {
        result(0) = betaDotN*psiVal*Leg2D(xMinus, yMinus, m, UMinus);
        return result;
    } else
    {
        // negative iPlus indicates exterior boundary
        if (iPlus < 0)
        {
            result(0) = 0;
            //dynamic_cast<PeriodicMesh &>(msh).getPeriodicCoordinates(iMinus, iPlus, x, y, xPlus, yPlus);
            //result(0) = betaDotN*psiVal*Leg2D(xPlus, yPlus, m, UPlus);
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
    
    result(0) = Leg2D(xx, yy, m, phi)*(beta_x(x, y)*Leg2D(xx, yy, m, psi_x) 
                                     + beta_y(x, y)*Leg2D(xx, yy, m, psi_y));
    return result;
}

mat Advection::computeBoundaryJacobian(double x, double y)
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
        
        result(0) = betaDotN*Leg2D(xPsi, yPsi, m, psi)*Leg2D(xPhi, yPhi, m, phi);
        
    } else
    {
        result(0) = 0;
    }
    
    return result;
}
