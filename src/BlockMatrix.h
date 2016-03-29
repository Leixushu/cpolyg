#pragma once

#include <armadillo>
#include <vector>

/// Abstract preconditioner class
struct Preconditioner;

/// Sparse matrix stored in block compressed sparse row format
struct BlockMatrix
{
    /// number of (block) rows
    int n_rows;
    /// size of each block
    int bl;
    /// total number of blocks
    int nb;
    
    /// list of all the blocks
    std::vector<arma::mat> blocks;
    /// list of the column indices for each block
    std::vector<int> colIndices;
    /// list consisting of the first index into blocks for each row
    std::vector<int> rowBlock;
    
    /// compute Ab = x given b, overwrite input array
    void matvec(double *b);
    /// compute Ab = x given b, return result
    arma::vec matvec(const arma::vec &b);
    
    /// Solve the system Ax = b using GMRES, initial guess provided by x.
    /** \param[in    ] b       right-hand side vector
        \param[in,out] x       on input initial guess vector, on output solution vector
        \param[in    ] m       restart parameter
        \param[in,out] tol     on input desired tolerance, on output norm of the residual
        \param[in,out] maxIt   on input maximum number of iterations, 
                               on output number of iterations used
        \param[in    ] pc      Preconditioner object to use
    */
    void gmres(arma::vec &b, arma::vec &x, int m, double &tol, int &maxIt, 
               Preconditioner &pc);
    
    /// Solve the system Ax = b using Jacobi iterations, return the result.
    /** \param[in    ] b        right-hand side vector
        \param[in,out] x        on input initial guess vector, on output solution vector
        \param[in,out] tol      on input desired tolerance, on output norm of the residual
        \param[in,out] maxIt    on input maximum number of iterations,
                                on output iterations used
        \param[in ] pc          BlockJacobi preconditioner object of the matrix
    */
    void jacobi(arma::vec &b, arma::vec &x, double &tol, int &maxIt, Preconditioner &pc);
    
    /// Output matrix in sparse (coordinate) format, suitable for plotting with gnuplot
    void spy(std::string filename);
    
    /// Scalar multiplication
    BlockMatrix& operator *=(double scale);
    BlockMatrix operator *(double scale);
    
    /// return BlockMatrix consisting of the diagonal blocks of M, block size a_b
    static BlockMatrix diag(arma::mat M, int a_b);
    /// return BlockMatrix consisting of the diagonal blocks of M
    static BlockMatrix diag(BlockMatrix &M);
    /// return BlockMatrix consisting of the off-diagonal blocks of M
    static BlockMatrix offDiag(BlockMatrix &M);
    
    static BlockMatrix blockBlockMatrix(arma::field<BlockMatrix> &b);
    
    /// default constructor
    BlockMatrix() {};
    
    BlockMatrix(int a_bl, int a_n);
};
