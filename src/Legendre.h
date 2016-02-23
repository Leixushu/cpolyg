#pragma once

#include <cstddef>
#include <gsl/gsl_sf_legendre.h>
#include <armadillo>

/*
Legendre Polynomial Routines
Will Pazner
November 18, 2015

---

Routines for compute values of a given Legendre series, and for computing the Legendre 
series coefficients of the x and y derivatives, given the coefficients of a series.
These routines can be compiled into a shared library, and then called from Python
using ctypes. The function leg2d makes use of the GNU Scientific Library function 
gsl_sf_legendre_Pl.

Example:
Degree 3 polynomials have coefficients c_{ij}, where i+j <= 3. We store these coefficients
as a one-dimensional array in the following format:

c0  c1  c2  c3
c4  c5  c6
c7  c8
c9

where the coefficients correspond to the basis polynomials:

1        L_1(y)         L_2(y)         L_3(y)
L_1(x)   L_1(x)L_1(y)   L_1(x)L_2(y) 
L_2(x)   L_2(x)L_1(y)
L_3(x)
*/

/*
Leg2D: compute the value of the bivariate Legendre series
inputs:
    double x   : x coordinate
    double y   : y coordinate
    int m      : degree of polynomial + 1
    double * c : coefficients of the Legendre series (length m(m+1)/2)
outputs:
    returns value:
        \sum_{i = 0}^{m-1} \sum_{i=0}^{m-1-j} c_{ij} L_i(x) L_j(y)

*/

// evaluate Legendre polynomial at a point (inlined for speed)
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
