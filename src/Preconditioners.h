#pragma once

#include "BlockMatrix.h"
#include <armadillo>
#include <vector>

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
    BlockILU0(BlockMatrix &M);
};
