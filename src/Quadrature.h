#pragma once

#include <armadillo>
#include "PolyMesh.h"

// Provide quadrature points and weights for various degree polynomials
// For integrals over lines and triangles
struct Quadrature
{
    static arma::mat lin2;
    static arma::mat lin4;
    static arma::mat lin6;
    static arma::mat lin10;
    
    static arma::mat tri2;
    static arma::mat tri4;
    static arma::mat tri6;
    static arma::mat tri10;
    
    static double polygonIntegral(const PolyMesh &msh, FnFunctor &cb, int pi, int deg);
    static arma::mat polygonIntegral(const PolyMesh &msh, VecFunctor &cb, int pi, int deg);
    
    static double lineIntegral(const PolyMesh &msh, FnFunctor &cb, int a, int b, int deg);
    static arma::mat lineIntegral(const PolyMesh &msh, VecFunctor &cb, int a, int b, int deg);
};
