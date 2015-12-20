#include "EulerVortex.h"
#include <cmath>

using namespace arma;
using namespace std;

double EulerVortex::ExactSolution::f(double x, double y) const
{
    return (1 - (pow((x - x0) - uBar*t, 2)) - pow((y - y0) - vBar*t,2))/(r_c*r_c);
}

double EulerVortex::ExactSolution::u(double x, double y) const
{
	return uInf*(cos(theta) - epsilon*((y - y0) - vBar*t)/(2*M_PI*r_c)*exp(f(x, y)*0.5));
}

double EulerVortex::ExactSolution::v(double x, double y) const
{
	return uInf*(sin(theta) + epsilon*((x - x0) - uBar*t)/(2*M_PI*r_c)*exp(f(x, y)*0.5));
}

double EulerVortex::ExactSolution::rho(double x, double y) const
{
	return rhoInf*pow(1 - epsilon*epsilon*(gamma - 1)*(MInf*MInf/(8*M_PI*M_PI))
                        *exp(f(x,y)), 1/(gamma-1));
}

double EulerVortex::ExactSolution::p(double x, double y) const
{
	return pInf*pow(1 - epsilon*epsilon*(gamma - 1)*(MInf*MInf/(8.0*M_PI*M_PI))
                        *exp(f(x,y)), (1/(gamma-1)));
}

double EulerVortex::ExactSolution::rhoE(double x, double y) const
{
	return p(x, y)/(gamma - 1) + 0.5*rho(x, y)*(pow(u(x,y),2) + pow(v(x,y),2));
}

vec EulerVortex::ExactSolution::operator()(double x, double y) const
{
    vec result(kEulerComponents);
    
    result[0] = rho(x, y);
    result[1] = rho(x, y)*u(x,y);
    result[2] = rho(x, y)*v(x,y);
    result[3] = rhoE(x, y);
    
    return result;
}

EulerVortex::EulerVortex(PolyMesh &m, double g) : Euler(m, g)
{
    exact.nc = kEulerComponents;
    
    exact.gamma = gamma;
    exact.x0 = 5;
    exact.y0 = 5;
    exact.theta = atan2(1, 2);
    exact.epsilon = 0.3;
    exact.r_c = 1.5;
    exact.MInf = 0.5;
    exact.uInf = 1.0;
    exact.rhoInf = 1.0;
    
    exact.uBar = exact.uInf*cos(exact.theta);
    exact.vBar = exact.uInf*sin(exact.theta);
    exact.EInf = (1.0/(gamma*exact.MInf*exact.MInf*(gamma-1)) + 0.5)*exact.uInf*exact.uInf;
    exact.pInf = (gamma-1.0)*exact.rhoInf*(exact.EInf - 0.5*exact.uInf*exact.uInf);
}

MeshFn EulerVortex::exactSolution(const double t, const int deg)
{
    MeshFn result(msh, deg, kEulerComponents);
    
    exact.t = t;
    
    result.interp(exact);
    
    return result;
}
