#include "KelvinHelmholtz.h"
#include <cmath>

using namespace arma;
using namespace std;

double KelvinHelmholtz::KHSolution::u(double x, double y) const
{
	return rho(x,y) - 1;
}

double KelvinHelmholtz::KHSolution::v(double x, double y) const
{
	return 0;
}

double KelvinHelmholtz::KHSolution::rho(double x, double y) const
{
	if ( fabs(y - 0.5) < (0.15 + sin(2*M_PI*x)/200.0) )
	{
	    return 2;
	} else
	{
	    return 1;
	}
}

double KelvinHelmholtz::KHSolution::p(double x, double y) const
{
	return 3;
}

double KelvinHelmholtz::KHSolution::rhoE(double x, double y) const
{
	return p(x, y)/(gamma - 1) + 0.5*rho(x, y)*(pow(u(x,y),2) + pow(v(x,y),2));
}

mat KelvinHelmholtz::KHSolution::operator()(double x, double y) const
{
    vec result(kEulerComponents);
    
    result[0] = rho(x, y);
    result[1] = rho(x, y)*u(x,y);
    result[2] = rho(x, y)*v(x,y);
    result[3] = rhoE(x, y);
    
    return result;
}

KelvinHelmholtz::KelvinHelmholtz(PeriodicMesh &m, double g) : Euler(m, g)
{
    KHSolution *kh = new KHSolution;
    
    kh->gamma = gamma;
    kh->t = 0;
    
    exact = kh;
}

MeshFn KelvinHelmholtz::exactSolution(const double t, const int deg)
{
    MeshFn result(msh, deg, kEulerComponents);
    
    exact->t = t;
    result.interp(*exact);
    
    return result;
}
