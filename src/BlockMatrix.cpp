#include "BlockMatrix.h"
#include "Preconditioners.h"
#include "blas/blas.h"
#include "Timer/CH_Timer.H"
#include <cassert>

using namespace std;
using namespace arma;

vec BlockMatrix::matvec(const vec &b) const
{
    vec result = zeros<vec>(b.n_rows);
    int i, j, k;
    
    for (i = 0; i < n_rows; i++)
    {
        for (k = rowBlock[i]; k < rowBlock[i+1]; k++)
        {
            j = colIndices[k];
            
            result.subvec(i*bl, (i+1)*bl - 1) += blocks[k]*b.subvec(j*bl, (j+1)*bl - 1);
        }
    }
    
    return result;
}

void BlockMatrix::matvec(double *b) const
{
    double *bCopy = new double[n_rows*bl];
    int i, j, k;
    
    copy(b, b+n_rows*bl, bCopy);
    for (i = 0; i < n_rows*bl; i++)
        b[i] = 0;
    
    for (i = 0; i < n_rows; i++)
    {
        for (k = rowBlock[i]; k < rowBlock[i+1]; k++)
        {
            j = colIndices[k];

            cdgemv('N', bl, bl, 1, const_cast<double *>(blocks[k].memptr()), bl,
                   bCopy + j * bl, 1.0, 1.0, b + i * bl, 1.0);
        }
    }
    
    delete[] bCopy;
}

void BlockMatrix::jacobi(arma::vec &b, arma::vec &x, double &tol, int &maxIt, Preconditioner &pc)
{
    vec Dinvb = b;
    vec ur;
    vec DinvAur;
    vec r;
    int m;
    double omega = 1.0;
    
    pc.solve(Dinvb.memptr());
    
    for (m = 0; m < maxIt; m++)
    {
        DinvAur = x;
        matvec(DinvAur.memptr());
        
        r = DinvAur - b;
        if (norm(r) < tol) break;
        
        pc.solve(DinvAur.memptr());
        
        x = omega*Dinvb + x - omega*DinvAur;
    }
    
    tol = norm(r);
    maxIt = m;
}

void BlockMatrix::gmres(vec &bvec, vec &xvec, int m, double &tol, int &maxit, Preconditioner &pc)
{
    CH_TIMERS("gmres");
    int n = n_rows*bl;
    double *b  = bvec.memptr();
    double *x  = xvec.memptr();
    double *v  = new double[n*(m+1)];
    double *h  = new double[m*(m+1)/2];
    double *r  = new double[n];
    double *y  = new double[m+1];
    double *c  = new double[m];
    double *s  = new double[m];

    double beta, hnew, rd, dd, nrm2b;
    int i=0, j, uki, u0i;

    copy(b,b+n,r); 
    pc.solve(r);
    nrm2b = cdnrm2(n,r,1);
    
    if (maxit%m>0)
    maxit=(maxit/m+1)*m;

    for (j=0; j<maxit/m; j++)
    {
        copy(x,x+n,r);
        matvec(r);
        cdaxpy(n,-1.,b,1,r,1);
        pc.solve(r);
        beta=cdnrm2(n,r,1);
        cdcopy(n,r,1,v,1);
        cdscal(n,1./beta,v,1);

        y[0]=beta;
        uki=0;
        for (i=0; i<m; i++)
        {
            if (fabs(y[i])<tol*nrm2b) break;
            double *vi=v+n*i;
            double *vi1=v+n*(i+1);
            u0i=uki;
            copy(vi,vi1,vi1);
            matvec(vi1);
            pc.solve(vi1);
            
            cdgemv('T',n,i+1,1.,v,n,vi1,1,0.,h+u0i,1);
            cdgemv('N',n,i+1,-1.,v,n,h+u0i,1,1.,vi1,1);
            
            hnew=cdnrm2(n,vi1,1);
            cdscal(n,1./hnew,vi1,1);
            for ( int k=0; k<i; ++k )
            {
                double tmp = c[k]*h[uki]-s[k]*h[uki+1];
                h[uki+1]   = s[k]*h[uki]+c[k]*h[uki+1];
                h[uki]     = tmp;
                ++uki;
            }

            rd     = h[uki];
            dd     = sqrt(rd*rd+hnew*hnew);
            c[i]   = rd/dd;
            s[i]   = -hnew/dd;
            h[uki] = dd;
            ++uki;

            y[i+1] = s[i]*y[i];
            y[i]   = c[i]*y[i];
        }

        cdtpsv('U','N','N',i,h,y,1);
        cdgemv('N',n,i,-1.,v,n,y,1,1.,x,1);

        if ( fabs(y[i])<tol*nrm2b ) break;
    }
    
    maxit = m*j+i;
    tol=std::abs(y[i])/nrm2b;

    delete[] v;
    delete[] h;
    delete[] r;
    delete[] y;
    delete[] c;
    delete[] s;
}

// construct a block-diagonal matrix from the diagonal blocks of M of size a_b
BlockMatrix BlockMatrix::diag(arma::mat M, int a_b)
{
    BlockMatrix result;
    mat block;
    int i;
    int b = a_b;
    // number of rows of blocks is number of rows in matrix divided by block size
    int n = M.n_rows/b;
    
    result = BlockMatrix(b, n);

    result.nb = n;
    result.colIndices.resize(n);
    result.rowBlock.resize(n+1);
    
    for (i = 0; i < n; i++)
    {
        block = M.submat(i*b, i*b, (i+1)*b - 1, (i+1)*b - 1);
        result.blocks.push_back(block);
        
        result.colIndices[i] = i;
        result.rowBlock[i] = i;
    }
    
    result.rowBlock[n] = n;
    
    return result;
}

// extract the diagonal blocks from a BlockMatrix
BlockMatrix BlockMatrix::diag(BlockMatrix &M)
{
    BlockMatrix result;
    int i, k;
    int n = M.n_rows;
    
    result = BlockMatrix(M.bl, n);
    result.colIndices.resize(n);
    
    for (i = 0; i < n; i++)
    {
        result.rowBlock[i] = result.nb;
        // search the blocks in the row for the column on the diagonal
        for (k = M.rowBlock[i]; k < M.rowBlock[i+1]; k++)
        {
            // if we're on the diagonal, add this block
            if (M.colIndices[k] == i)
            {
                result.blocks.push_back(M.blocks[k]);
                result.colIndices[result.nb] = i;
                result.nb++;
                break;
            }
        }
    }
    
    // make sure to end the rowBlock vector
    result.rowBlock[n] = result.nb;
    
    return result;
}

BlockMatrix BlockMatrix::blockBlockMatrix(arma::field<BlockMatrix> &b)
{
    assert(b.n_rows == b.n_cols);
    int nBig = b.n_rows;
    int nSmall = b(0,0).n_rows;
    int bl = b(0,0).bl;
    int k;
    BlockMatrix result(bl, nBig*nSmall);
    
    for (size_t i = 0; i < nBig; i++)
    {
        for (size_t iSmall = 0; iSmall < nSmall; iSmall++)
        {
            result.rowBlock[i*nSmall + iSmall] = result.nb;
            for (size_t j = 0; j < nBig; j++)
            {
                for (k = b(i,j).rowBlock[iSmall]; 
                     k < b(i,j).rowBlock[iSmall+1]; k++)
                {
                    result.blocks.push_back(b(i,j).blocks[k]);
                    result.colIndices.push_back(j*nSmall + b(i,j).colIndices[k]);
                    result.nb++;
                }
            }
        }
    }
    
    result.rowBlock[result.n_rows] = result.nb;
    
    return result;
}

// extract the off-diagonal blocks from a BlockMatrix
BlockMatrix BlockMatrix::offDiag(BlockMatrix &M)
{
    BlockMatrix result;
    int i, k, j;
    int n = M.n_rows;
    
    // use the same block size, and same number of rows
    result = BlockMatrix(M.bl, n);
    
    for (i = 0; i < n; i++)
    {
        result.rowBlock[i] = result.nb;
        
        // search the blocks in this row for off-diagonal columns
        for (k = M.rowBlock[i]; k < M.rowBlock[i+1]; k++)
        {
            j = M.colIndices[k];
            
            // if we're not on the diagonal
            if (i != j)
            {
                // add the block
                result.blocks.push_back(M.blocks[k]);
                result.colIndices.push_back(j);
                result.nb++;
            }
        }
    }
    
    // make sure to end the rowBlock vector
    result.rowBlock[n] = result.nb;
    
    return result;
}


BlockMatrix& BlockMatrix::operator *=(double scale)
{
    int i;
    
    // multiply each block by a scalar
    for (i = 0; i < nb; i++)
    {
        blocks[i] *= scale;
    }
    
    return *this;
}

BlockMatrix BlockMatrix::operator *(double scale)
{
    BlockMatrix result = *this;
    result *= scale;
    return result;
}

// write the matrix to a file in a format readable by gnuplot
void BlockMatrix::spy(std::string filename)
{
    int i, j, k, ii, jj;
    ofstream file;
    file.open(filename);
    
    // loop over the rows
    for (i = 0; i < n_rows; i++)
    {
        // loop over the blocks in each row
        for (k = rowBlock[i]; k < rowBlock[i+1]; k++)
        {
            j = colIndices[k];
            
            // loop over the entries of each block
            for (ii = 0; ii < bl; ii++)
            {
                for (jj = 0; jj < bl; jj++)
                {
                    // write only the nonzero entries
                    if (fabs(blocks[k](ii, jj)) > kEPS)
                    {
                        file << i*bl + ii << "\t" << j*bl + jj << "\t"
                             << blocks[k](ii, jj) << endl;
                    }
                }
            }
        }
    }
    
    file.close();
}

BlockMatrix::BlockMatrix(int a_bl, int a_n)
{
    bl = a_bl;
    n_rows = a_n;
    nb = 0;
    
    blocks.clear();
    colIndices.clear();
    rowBlock = vector<int>(a_n + 1, 0);
}
