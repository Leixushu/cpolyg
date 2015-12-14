#include <gsl/gsl_sf_legendre.h>

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
double Leg2D(double x, double y, int m, double * c)
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

/*
LegDerX: compute the x derivative of the given bivariate Legendre series
inputs:
    int m       : degree of polynomial + c
    double * c  : Legendre coefficients (length m(m+1)/2)
outputs:
    double * cp : Legendre coefficients of the x derivative (must be already allocated)
*/
void LegDerX(int m, double * c, double * cp)
{
    int i, j, k;
    
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
}

/*
LegDerY: compute the y derivative of the given bivariate Legendre series
inputs:
    int m       : degree of polynomial + c
    double * c  : Legendre coefficients (length m(m+1)/2)
outputs:
    double * cp : Legendre coefficients of the y derivative (must be already allocated)
*/
void LegDerY(int m, double * c, double * cp)
{
    int i, j, k;
    
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
}
