#include "Preconditioners.h"
#include "blas/blas.h"

using namespace std;
using namespace arma;

void BlockJacobi::solve(double *x)
{
    int i;
    double *A = new double[bl*bl];
    int *ipvt = new int[bl];
    
    // do a block diagonal solve, overwrite the vector x
    for (i = 0; i < n; i++)
    {
        copy(blocks[i].memptr(), blocks[i].memptr() + bl*bl, A);
        // LAPACK LU factorization
        cdgetrf(bl, bl, A, bl, ipvt);
        // LU dense solve
        cdgetrs('N', bl, 1, A, bl, ipvt, x + i*bl, bl);
    }
    
    delete[] A;
    delete[] ipvt;
}

BlockJacobi::BlockJacobi(BlockMatrix &M)
{
    int i, j, k;
    blocks.resize(M.n_rows, zeros<mat>(M.bl, M.bl));
    
    n = M.n_rows;
    bl = M.bl;
    
    // extract the diagonal blocks from M
    for (i = 0; i < n; i++)
    {
        for (k = M.rowBlock[i]; k < M.rowBlock[i+1]; k++)
        {
            j = M.colIndices[k];
            if (j == i)
            {
                blocks[i] = M.blocks[k];
            }
        }
    }
}
