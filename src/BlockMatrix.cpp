#include "BlockMatrix.h"
#include "fortwrap.h"

using namespace std;
using namespace arma;

vec BlockMatrix::matvec(const vec &b)
{
    vec result = zeros<vec>(b.n_rows);
    int i, j, k;
    
    for (i = 0; i < nb; i++)
    {
        for (k = rowBlock[i]; k < rowBlock[i+1]; k++)
        {
            j = colIndices[k];
            
            result.subvec(i*bl, (i+1)*bl - 1) += blocks[k]*b.subvec(j*bl, (j+1)*bl - 1);
        }
    }
    
    return result;
}

void BlockMatrix::matvec(double *b)
{
    double *bCopy = new double[nb*bl];
    int i, j, k;
    
    copy(b, b+nb*bl, bCopy);
    for (i = 0; i < nb*bl; i++)
        b[i] = 0;
    
    for (i = 0; i < nb; i++)
    {
        for (k = rowBlock[i]; k < rowBlock[i+1]; k++)
        {
            j = colIndices[k];
            
            cdgemv('N', bl, bl, 1, blocks[k].memptr(), bl, bCopy + j*bl,
                    1.0, 1.0, b + j*bl, 1.0);
        }
    }
    
    delete[] bCopy;
}

void BlockMatrix::gmres(vec &bvec, vec &xvec, int m, double &tol, int &maxit, Preconditioner &pc)
{
    int n = nb*bl;
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
            // psolve(vi1,fpar);
            
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

BlockMatrix BlockMatrix::diag(arma::mat M, int a_b)
{
    BlockMatrix result;
    mat block;
    int i;
    int b = a_b;
    int n = M.n_rows/b;
    
    result.bl = a_b;
    result.nb = n;
    
    result.blocks.clear();
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
