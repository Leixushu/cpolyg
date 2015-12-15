#pragma once

#include <armadillo>
#include "PolyMesh.h"

typedef double (*FnCallback)(double x, double y);

struct PolyFn
{
    arma::vec a;
    
    double eval(double x, double y, PolyMesh &msh, int p);
};

struct MeshFn
{
    PolyMesh &msh;
    int nc;
    int deg;
    arma::cube a;
    
    MeshFn(PolyMesh &a_msh, int a_deg, int a_nc = 1);
    MeshFn(PolyMesh &a_msh, FnCallback cb, int a_deg);
    
    void interp(FnCallback cb, int c);
    
    double eval(double x, double y, PolyMesh &msh, int p);
    PolyFn getPolyFn(int p, int c);
};