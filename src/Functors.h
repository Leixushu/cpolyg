#pragma once

#include <armadillo>
#include <iostream>

typedef double (*FnCallback)(double x, double y);
typedef arma::vec (*VecCallback)(double x, double y);

/// Abstract class for a functor that takes (x,y) coordinates and return a function value
struct FnFunctor
{
    virtual double operator()(double x, double y) const = 0;
    
    FnFunctor() {};
    //FnFunctor(const FnFunctor &fn) { };
    
    virtual ~FnFunctor() {};
};

/// Functor that wraps a basic callback function
struct FnCallbackFunctor : FnFunctor
{   
    FnCallback cb;
    
    FnCallbackFunctor(FnCallback fn) : cb(fn) { }
    
    double operator()(double x, double y) const
    {
        return cb(x, y);
    }
};

/// Abstract class for a functor that takes (x,y) coordinates and returns a vector
struct VecFunctor
{
    int n_rows;
    int n_cols;
    
    VecFunctor() : n_rows(1), n_cols(1) { };
    
    virtual arma::mat operator()(double x, double y) const = 0;
    virtual ~VecFunctor() {};
};

struct VecCallbackFunctor : VecFunctor
{   
    VecCallback cb;
    
    VecCallbackFunctor(VecCallback fn, int a_nc) : cb(fn)
    {
        n_rows = a_nc;
    }
    
    arma::mat operator()(double x, double y) const
    {
        return cb(x, y);
    }
};

struct VecFnFunctor : VecFunctor
{
    const FnFunctor &fncb;
    VecFnFunctor(const FnFunctor &fn) : fncb(fn) { };
    
    arma::mat operator()(double x, double y) const
    {
        return arma::vec::fixed<1>({fncb(x, y)});
    }
};
