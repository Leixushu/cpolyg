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
    MeshFn(const MeshFn &fn);
    
    void interp(const FnFunctor &cb, int component);
    void interp(const VecFunctor &cb, int component = 0);
    double eval(double x, double y, int p, int c = 0) const;
    
    void gnuplot(std::string filename) const;
    
    MeshFn& operator+=(const MeshFn fn2)
    {
        a += fn2.a;
        return (*this);
    }
    MeshFn operator+(const MeshFn &fn2) const;
    MeshFn& operator=(const MeshFn &fn);
    MeshFn operator*(const double scale) const;
    MeshFn operator/(const double scale) const
    {
        return (*this)*(1.0/scale);
    }
};

inline MeshFn operator*(const double scale, const MeshFn& fn)
{
    return fn*scale;
}
