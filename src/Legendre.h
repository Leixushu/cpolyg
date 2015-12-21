#pragma once

#include <gsl/gsl_sf_legendre.h>
#include <armadillo>


inline double Leg2D(double x, double y, int m, arma::vec c)
{
    int i, j, k;
    double val = 0;
    
    k = 0;
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m-i; j++)
        {
            val += c[k]*gsl_sf_legendre_Pl(i, x)*gsl_sf_legendre_Pl(j, y);
            k++;
        }
    }
    
    return val;
}

arma::vec LegDerX(int m, arma::vec &c);
arma::vec LegDerY(int m, arma::vec &c);
