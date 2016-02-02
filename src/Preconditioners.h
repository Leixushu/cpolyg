#pragma once

#include <armadillo>
#include <vector>
#include "BlockMatrix.h"
#include "Jacobian.h"

struct NoPreconditioner : Preconditioner
{
    void solve(double *x) { };
};

struct BlockJacobi : Preconditioner
{
    int n;
    int bl;
    std::vector<arma::mat> blocks;
    
    void solve(double *x);
    
    BlockJacobi(BlockMatrix &M);
};

struct BlockILU0 : Preconditioner
{
    int n;
    int bl;
    
    BlockMatrix DD, L, DU;
    BlockILU0(Jacobian &A);
    
    void solve(double *x);
};
