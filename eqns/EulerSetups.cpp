#include "EulerSetups.h"
#include <cmath>

using namespace arma;
using namespace std;

double VortexSolution::f(double x, double y) const
{
    return (1 - (pow((x - x0) - uBar*t, 2)) - pow((y - y0) - vBar*t,2))/(r_c*r_c);
}

double VortexSolution::u(double x, double y) const
{
	return uInf*(cos(theta) - epsilon*((y - y0) - vBar*t)/(2*M_PI*r_c)*exp(f(x, y)*0.5));
}

double VortexSolution::v(double x, double y) const
{
	return uInf*(sin(theta) + epsilon*((x - x0) - uBar*t)/(2*M_PI*r_c)*exp(f(x, y)*0.5));
}

double VortexSolution::rho(double x, double y) const
{
	return rhoInf*pow(1 - epsilon*epsilon*(gamma - 1)*(MInf*MInf/(8*M_PI*M_PI))
                        *exp(f(x,y)), 1/(gamma-1));
}

double VortexSolution::p(double x, double y) const
{
	return pInf*pow(1 - epsilon*epsilon*(gamma - 1)*(MInf*MInf/(8.0*M_PI*M_PI))
                        *exp(f(x,y)), (1/(gamma-1)));
}

double VortexSolution::rhoE(double x, double y) const
{
	return p(x, y)/(gamma - 1) + 0.5*rho(x, y)*(pow(u(x,y),2) + pow(v(x,y),2));
}

mat VortexSolution::operator()(double x, double y) const
{
    vec::fixed<kEulerComponents> result;
    
    result[0] = rho(x, y);
    result[1] = rho(x, y)*u(x,y);
    result[2] = rho(x, y)*v(x,y);
    result[3] = rhoE(x, y);
    
    return result;
}

VortexSolution::VortexSolution(double a_gamma) : VecFunctor(kEulerComponents)
{
    gamma = a_gamma;
    x0 = 5;
    y0 = 5;
    theta = atan2(1, 2);
    epsilon = 0.3;
    r_c = 1.5;
    MInf = 0.5;
    uInf = 1.0;
    rhoInf = 1.0;
    
    uBar = uInf*cos(theta);
    vBar = uInf*sin(theta);
    
    pInf = 1;
}

MeshFn VortexSolution::interpolated(const PolyMesh &msh, const int deg)
{
    MeshFn result(msh, deg, kEulerComponents);
    result.interp(*this);
    return result;
}

double KelvinHelmholtz::phi(double y) const
{
    return y/10.0 - (0.4)*(erf(26.66*(y-4))-1);
}

double KelvinHelmholtz::psi(double x) const
{
    return 1 + cos(M_PI*x/20.0);
}

double KelvinHelmholtz::u(double x, double y) const
{
	return rho(x,y) - 1;
	//return sqrt(gamma)*psi(x);
}

double KelvinHelmholtz::v(double x, double y) const
{
	return 0;
}

double KelvinHelmholtz::rho(double x, double y) const
{
	if ( fabs(y - 0.5) < (0.15 + sin(2*M_PI*x)/200.0) )
	{
	    return 2;
	} else
	{
	    return 1;
	}
    //return 1 + 0.01*psi(x) + phi(y);
}

double KelvinHelmholtz::p(double x, double y) const
{
	return 3;
	//return (1 + gamma/20.0*psi(x))*400.0;
}

double KelvinHelmholtz::rhoE(double x, double y) const
{
	return p(x, y)/(gamma - 1) + 0.5*rho(x, y)*(pow(u(x,y),2) + pow(v(x,y),2));
}

mat KelvinHelmholtz::operator()(double x, double y) const
{
    vec::fixed<kEulerComponents> result;
    
    result[0] = rho(x, y);
    result[1] = rho(x, y)*u(x,y);
    result[2] = rho(x, y)*v(x,y);
    result[3] = rhoE(x, y);
    
    return result;
}

MeshFn KelvinHelmholtz::interpolated(const PolyMesh &msh, const int deg)
{
    MeshFn result(msh, deg, kEulerComponents);
    result.interp(*this);
    return result;
}
