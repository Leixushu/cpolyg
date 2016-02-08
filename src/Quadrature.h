#pragma once

#include <armadillo>

// Provide quadrature points and weights for various degree polynomials
// For integrals over lines and triangles
struct Quadratures
{
    static arma::mat lin2;
    static arma::mat lin4;
    static arma::mat lin6;
    static arma::mat lin10;
    
    static arma::mat tri2;
    static arma::mat tri4;
    static arma::mat tri6;
    static arma::mat tri10;
};
