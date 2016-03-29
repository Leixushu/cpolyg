#include "Preconditioners.h"
#include "blas/blas.h"
#include "Timer/CH_Timer.H"

using namespace std;
using namespace arma;

void BlockJacobi::solve(double *x)
{
    int i;
    // do a block diagonal solve, overwrite the vector x
    for (i = 0; i < n; i++)
    {
        // Solve using the LU factorization of the diagonal blocks
        cdgetrs('N', bl, 1, factoredBlocks(i).memptr(), bl, pivots(i).memptr(), 
                x + i*bl, bl);
    }
}

BlockJacobi::BlockJacobi(BlockMatrix &M)
{
    int i, j, k;
    
    n = M.n_rows;
    bl = M.bl;
    
    factoredBlocks.set_size(n);
    pivots.set_size(n);
    
    // extract the diagonal blocks from M
    for (i = 0; i < n; i++)
    {
        for (k = M.rowBlock[i]; k < M.rowBlock[i+1]; k++)
        {
            j = M.colIndices[k];
            if (j == i)
            {
                // allocate memory for the pivot indices
                pivots(i).set_size(bl);
                factoredBlocks(i) = M.blocks[k];
                // LAPACK LU factorization
                cdgetrf(bl, bl, factoredBlocks(i).memptr(), bl, pivots(i).memptr());
                break;
            }
        }
    }
}

BlockILU0::BlockILU0(BlockMatrix &A)
{
    CH_TIMERS("Compute ILU(0)");
    int i, j, k, k2;
    mat D, Id;
    BlockMatrix AD, AO;
    
    n = A.n_rows;
    bl = A.bl;
    
    // initialize the block matrix consisting of diagonal blocks of A
    AD = BlockMatrix::diag(A);
    AO = BlockMatrix::offDiag(A);
    
    DD = BlockMatrix(bl, n);
    DD.colIndices.resize(n);
    
    // loop over all rows of A
    for (i = 0; i < A.n_rows; i++)
    {
        // extract the diagonal block and invert it
        D = inv(AD.blocks[AD.rowBlock[i]]);
        DD.blocks.push_back(D);
        DD.colIndices[i] = i;
        DD.rowBlock[i] = i;
        
        // loop over all non-diagonal blocks
        for (k = AO.rowBlock[i]; k < AO.rowBlock[i+1]; k++)
        {
             // get the column index
            j = AO.colIndices[k];
            
            if (j > i)
            {
                for (k2 = AO.rowBlock[j]; k2 < AO.rowBlock[j+1]; k2++)
                {
                    if (AO.colIndices[k2] == i)
                    {
                        // do increment operators here -- probably faster than assignment
                        AO.blocks[k2] = AO.blocks[k2]*D;
                        AD.blocks[j] = AD.blocks[j] - AO.blocks[k2]*AO.blocks[k];
                        AO.blocks[k] = D*AO.blocks[k];
                        break;
                    }
                }
            }
        }
    }
    DD.rowBlock[n] = n;
    
    // now get the upper and lower parts of AO
    L = BlockMatrix(bl, n);
    DU = BlockMatrix(bl, n);
    
    // loop over rows of AO
    for (i = 0; i < AO.n_rows; i++)
    {
        // normally L and DU would have identity matrices along the diagonal
        // but we use them implicitly in the backward and forward triangular solvers
        L.rowBlock[i] = L.nb;
        DU.rowBlock[i] = DU.nb;
        
        for (k = AO.rowBlock[i]; k < AO.rowBlock[i+1]; k++)
        {
            j = AO.colIndices[k];
            
            // if the block is to the left of the diagonal, add it to L
            // otherwise add it to DU
            if( j < i)
            {
                L.blocks.push_back(AO.blocks[k]);
                L.colIndices.push_back(j);
                L.nb++;
            } else if (j > i)
            {
                DU.blocks.push_back(AO.blocks[k]);
                DU.colIndices.push_back(j);
                DU.nb++;
            }
        }
    }
    
    L.rowBlock[n] = L.nb;
    DU.rowBlock[n] = DU.nb;
}

void BlockILU0::solve(double *x)
{
    int i, j, k;
    
    // do forward substitution to find L y = x (overwrite x with y)
    // L equal to the identity matrix + strictly lower triangular
    // so first block of y = first block of x
    for (i = 0; i < n; i++)
    {
        for (k = L.rowBlock[i]; k < L.rowBlock[i+1]; k++)
        {
            j = L.colIndices[k];
            // subtract A_j y from the right hand for every block left of the diagonal
            // i.e. y = x - sum_{j < i} A_j y_j
            // we do this by setting alpha = 1 and beta = 1 in the DGEMV call
            cdgemv('N', bl, bl, -1, L.blocks[k].memptr(), bl, x+bl*j, 1, 1, x+bl*i, 1);
        }
    }
    
    // multiply by D
    DD.matvec(x);
    
    // do backward substitution to find DU y = x (overwrite x with y)
    // DU is equal to identity matrix + strictly upper triangular (analogous to the
    // previous case, but loop over rows from last to first)
    for (i = n-1; i >= 0; i--)
    {
        for (k = DU.rowBlock[i]; k < DU.rowBlock[i+1]; k++)
        {
            j = DU.colIndices[k];
            cdgemv('N', bl, bl, -1, DU.blocks[k].memptr(), bl, x+bl*j, 1, 1, x+bl*i, 1);
        }
    }
}
