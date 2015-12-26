#pragma once

#include <armadillo>
#include <vector>

struct Preconditioner
{
    virtual void solve(double *x) = 0;
};

struct BlockMatrix
{
    // size of each block
    int bl;
    // number of blocks
    int nb;
    // number of (block) rows
    int n_rows;
    
    std::vector<arma::mat> blocks;
    std::vector<int> colIndices;
    std::vector<int> rowBlock;
    
    void matvec(double *b);
    arma::vec matvec(const arma::vec &b);
    
    void gmres(arma::vec &b, arma::vec &x, int m, double &tol, int &maxIt, 
               Preconditioner &pc);
    arma::vec jacobi(arma::vec &b, double &tol, int &maxIt, Preconditioner &pc);
    
    void spy(std::string filename);
    
    static BlockMatrix diag(arma::mat M, int a_b);
    
    BlockMatrix& operator *=(double scale);
    
    BlockMatrix() {};
};
