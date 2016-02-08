#pragma once

#include <armadillo>
#include "BlockMatrix.h"
#include "Jacobian.h"

// Provide various preconditioners for use in linear solvers

// Abstract class for preconditioners
struct Preconditioner
{
    // Perform solve Mx = b, overwriting input array
    virtual void solve(double *x) = 0;
};

// This 'preconditioner' does nothing
struct NoPreconditioner : Preconditioner
{
    void solve(double *x) { };
};

// Block Jacobi preconditioner: M = diagonal blocks of A
struct BlockJacobi : Preconditioner
{
    // number of block rows
    int n;
    // size of each block
    int bl;
    // store the LU factorization of each block (and pivot indices)
    arma::field<arma::mat> factoredBlocks;
    arma::field<arma::Col<int>> pivots;
    
    // do the block diagonal solve
    void solve(double *x);
    
    // extract the diagonal blocks and perform the factorization for each
    BlockJacobi(BlockMatrix &M);
};

// Block ILU(0) preconditioner: incomplete LU with same (block) sparsity as A
struct BlockILU0 : Preconditioner
{
    // number of block rows
    int n;
    // size of each block
    int bl;
    
    // store the diagonal, L, and U matrices in block compressed sparse row format
    BlockMatrix DD, L, DU;
    
    // do forward and then backward substitution: DU\(DD*(L\x))
    void solve(double *x);
    // perform the factorization
    BlockILU0(Jacobian &A);
};
