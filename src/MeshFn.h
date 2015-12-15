#pragma once

#include <armadillo>
#include "PolyMesh.h"

struct MeshFn
{
    PolyMesh &msh;
    int nc;
    int deg;
    arma::cube a;
    
    MeshFn(PolyMesh &a_msh, int a_deg, int a_nc = 1);
    MeshFn(PolyMesh &a_msh, FnCallback cb, int a_deg);
    
    void interp(FnCallback cb, int c);
    double eval(double x, double y, int p, int c = 0);
    
    void gnuplot(std::string filename, int c = 0);
    
    MeshFn operator+(MeshFn &fn2);
};
