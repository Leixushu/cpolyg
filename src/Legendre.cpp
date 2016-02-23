#include <stdlib.h>
#include <string.h>
#include "Legendre.h"

using namespace arma;

// inline double Leg2D(double x, double y, int m, vec c)
// {
//     int i, j, k;
//     double val = 0;
//     
//     k = 0;
//     
//     for (i = 0; i < m; i++)
//     {
//         for (j = 0; j < m-i; j++)
//         {
//             val += c[k]*gsl_sf_legendre_Pl(i, x)*gsl_sf_legendre_Pl(j, y);
//             k++;
//         }
//     }
//     
//     return val;
// }

/*
LegDerX: compute the x derivative of the given bivariate Legendre series
inputs:
    int m       : degree of polynomial + c
    double * c  : Legendre coefficients (length m(m+1)/2)
outputs:
    double * cp : Legendre coefficients of the x derivative (must be already allocated)
*/
vec LegDerX(int m, vec &a_c)
{
    int i, j, k;
    vec c = a_c;
    vec cp = zeros<vec>(m*(m+1)/2);
    
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m-i; j++)
        {
            cp[k] = 0;
            k++;
        }
    }
    
    // loop over the rows
    for (i = m - 1; i > 2; i--)
    {
        for (j = 0; j < m-i; j++)
        {
            cp[-(i-1)*(i-2*(1 + m))/2 + j] = (2*i - 1)*c[-i*(i - 1 - 2*m)/2 + j];
            c[-(i-2)*(i-3-2*m)/2 + j] += c[-i*(i - 1 - 2*m)/2 + j];
        }
    }
    
    if (m > 2)
    {
        for (j = 0; j < m - 2; j++)
        {
            cp[m+j] = 3*c[2*m - 1 + j];
        }
    }
    
    for (j = 0; j < m - 1; j++)
    {
        cp[j] = c[m+j];
    }
    
    return cp;
}

/*
LegDerY: compute the y derivative of the given bivariate Legendre series
inputs:
    int m       : degree of polynomial + c
    double * c  : Legendre coefficients (length m(m+1)/2)
outputs:
    double * cp : Legendre coefficients of the y derivative (must be already allocated)
*/
vec LegDerY(int m, vec &a_c)
{
    int i, j, k;
    vec c = a_c;
    vec cp = zeros<vec>(m*(m+1)/2);
    
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m-i; j++)
        {
            cp[k] = 0;
            k++;
        }
    }
    
    // loop over the columns
    for (j = m - 1; j > 2; j--)
    {
        for (i = 0; i < m-j; i++)
        {
            cp[-i*(i - 1 - 2*m)/2 + j - 1] = (2*j - 1)*c[-i*(i - 1 - 2*m)/2 + j];
            c[-i*(i - 1 - 2*m)/2 + j - 2] += c[-i*(i - 1 - 2*m)/2 + j];
        }
    }
    
    if (m > 2)
    {
        for (i = 0; i < m - 2; i++)
        {
            cp[-i*(i - 1 - 2*m)/2 + 1] = 3*c[-i*(i - 1 - 2*m)/2 + 2];
        }
    }
    
    for (i = 0; i < m - 1; i++)
    {
        cp[-i*(i - 1 - 2*m)/2] = c[-i*(i - 1 - 2*m)/2 + 1];
    }
    
    return cp;
}
