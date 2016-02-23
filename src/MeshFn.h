#pragma once

#include <armadillo>
#include "PolyMesh.h"

// Represent a function defined on a polygonal mesh, given by the coefficients 
// of its local basis expansion.
// The coefficients are represented as a cube of dimension (basisSize, nc, np)
// where nc is the number of components, and np is the number of polygons in the mesh.

struct MeshFn
{
    struct L2Difference : VecFunctor
    {
        const MeshFn &fn;
        const VecFunctor &exact;
        
        int i;
        
        L2Difference(const MeshFn &f, const VecFunctor &e)
        : VecFunctor(f.nc), fn(f), exact(e)
        { }
        
        arma::mat operator()(double x, double y) const;
    };
    
    struct L2Functor : VecFunctor
    {
        const MeshFn &fn;
        
        int i;
        
        L2Functor(const MeshFn &f) : VecFunctor(f.nc), fn(f) { }
        
        arma::mat operator()(double x, double y) const;
    };
    
    const PolyMesh &msh;
    int nc;
    int deg;
    arma::cube a;
    
    MeshFn(const PolyMesh &a_msh, int a_deg, int a_nc = 1);
    MeshFn(const PolyMesh &a_msh, FnCallback cb, int a_deg);
    MeshFn(const MeshFn &fn);
    
    void interp(const FnFunctor &cb, int component = 0);
    void interp(const VecFunctor &cb, int component = 0);
    double eval(double x, double y, int p, int c = 0) const;
    
    arma::vec L2Norm() const;
    
    double L2Error(const FnFunctor &exact) const;
    arma::vec L2Error(const VecFunctor &exact) const;
    
    void gnuplot(std::string filename) const;
    
    MeshFn& operator+=(const MeshFn &fn2)
    {
        a += fn2.a;
        return (*this);
    }
    
    MeshFn& operator-=(const MeshFn &fn2)
    {
        a -= fn2.a;
        return (*this);
    }
    
    MeshFn operator+(const MeshFn &fn2) const;
    MeshFn operator-(const MeshFn &fn2) const;
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
