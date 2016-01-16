#include "EulerVortex.h"
#include <cmath>

using namespace arma;
using namespace std;

double EulerVortex::VortexSolution::f(double x, double y) const
{
    return (1 - (pow((x - x0) - uBar*t, 2)) - pow((y - y0) - vBar*t,2))/(r_c*r_c);
}

double EulerVortex::VortexSolution::u(double x, double y) const
{
	return uInf*(cos(theta) - epsilon*((y - y0) - vBar*t)/(2*M_PI*r_c)*exp(f(x, y)*0.5));
}

double EulerVortex::VortexSolution::v(double x, double y) const
{
	return uInf*(sin(theta) + epsilon*((x - x0) - uBar*t)/(2*M_PI*r_c)*exp(f(x, y)*0.5));
}

double EulerVortex::VortexSolution::rho(double x, double y) const
{
	return rhoInf*pow(1 - epsilon*epsilon*(gamma - 1)*(MInf*MInf/(8*M_PI*M_PI))
                        *exp(f(x,y)), 1/(gamma-1));
}

double EulerVortex::VortexSolution::p(double x, double y) const
{
	return pInf*pow(1 - epsilon*epsilon*(gamma - 1)*(MInf*MInf/(8.0*M_PI*M_PI))
                        *exp(f(x,y)), (1/(gamma-1)));
}

double EulerVortex::VortexSolution::rhoE(double x, double y) const
{
	return p(x, y)/(gamma - 1) + 0.5*rho(x, y)*(pow(u(x,y),2) + pow(v(x,y),2));
}

mat EulerVortex::VortexSolution::operator()(double x, double y) const
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
    VortexSolution *vortex = new VortexSolution;
    
    vortex->gamma = gamma;
    vortex->x0 = 5;
    vortex->y0 = 5;
    vortex->theta = atan2(1, 2);
    vortex->epsilon = 0.3;
    vortex->r_c = 1.5;
    vortex->MInf = 0.5;
    vortex->uInf = 5.0;
    vortex->rhoInf = 1.0;
    
    vortex->uBar = vortex->uInf*cos(vortex->theta);
    vortex->vBar = vortex->uInf*sin(vortex->theta);
    
    vortex->pInf = 1;
    vortex->t = 0;
    
    exact = vortex;
    
    ((LaxFriedrichsFlux*)(boundaryTerm))->exact = exact;
}

MeshFn EulerVortex::exactSolution(const double t, const int deg)
{
    MeshFn result(msh, deg, kEulerComponents);
    
    exact->t = t;
    result.interp(*exact);
    
    return result;
}
