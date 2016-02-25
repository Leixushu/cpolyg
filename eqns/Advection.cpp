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

mat Advection::fluxFunction(const vec &vars, double x, double y)
{
    rowvec::fixed<2> result = {betaX(x, y)*vars(0), betaY(x, y)*vars(0)};
    return result;
}

vec Advection::numericalFluxFunction(const vec &varsMinus, const vec &varsPlus, 
    double x, double y, double nx, double ny)
{
    vec::fixed<1> result;
    double betaDotN = betaX(x, y)*nx + betaY(x, y)*ny;
    
    if (betaDotN > 0)
    {
        result(0) = varsMinus(0)*betaDotN;
    } else
    {
        result(0) = varsPlus(0)*betaDotN;
    }
    return result;
}

cube Advection::fluxJacobian(const vec &vars, double x, double y)
{
    cube::fixed<1,1,2> result;
    
    result(0,0,0) = betaX(x, y);
    result(0,0,1) = betaY(x, y);
    
    return result;
}

mat Advection::numericalFluxJacobian(const vec &varsMinus, const vec &varsPlus, double x, 
    double y, double nx, double ny, int sgn)
{
    vec::fixed<1> result;
    double betaDotN;
    
    betaDotN = nx*betaX(x, y) + ny*betaY(x, y);
    
    if (sgn*betaDotN > 0)
    {
        result = betaDotN;
    } else
    {
        result(0) = 0;
    }
    return result;
}
