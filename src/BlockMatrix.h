#pragma once

#include <armadillo>
#include <vector>

struct Preconditioner
{
    virtual void solve(double *x) = 0;
};

struct NoPreconditioner : Preconditioner
{
    void solve(double *x) { }
};

struct BlockMatrix
{
    // size of each block
    int bl;
    // number of blocks
    int nb;
    
    std::vector<arma::mat> blocks;
    std::vector<int> colIndices;
    std::vector<int> rowBlock;
    
    void matvec(double *b);
    arma::vec matvec(const arma::vec &b);
    void gmres(arma::vec &b, arma::vec &x, int m, double &tol, int &maxIt, 
               Preconditioner &pc);
    
    static BlockMatrix diag(arma::mat M, int a_b);
};
