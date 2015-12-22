#include "fortwrap.h"

/*
  General variable names:
 
  - lda, ldb, etc are leading dimensions of matrices A, B, etc
  - incx, incy, etc are increments of elements in x, y, etc
*/

//////////////////// DOUBLE PRECISION /////////////////////////////////

void cdlaswp(szint n, double *A, szint lda, szint k1, szint k2,
             szint *ipiv, szint incx)
  /* Perform row interchanges on A, from k1 to k2 in ipiv */
{
  DLASWP(&n,A,&lda,&k1,&k2,ipiv,&incx);
}

void cdtrsm(char side, char uplo, char transa, char diag, szint m, szint n,
            double alpha, double *A, szint lda, double *B, szint ldb)
  /* Solve AX=alpha*B or XA=alpha*B for X with A triangular */
{
  DTRSM(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}

void cdtpsv(char uplo, char trans, char diag, szint n, double *Ap, double *x, szint incx)
  /* Solve Ax=b, A is n-by-n triangular */ 
{
  DTPSV(&uplo,&trans,&diag,&n,Ap,x,&incx);
}

szint cdgesv(szint n, szint nrhs, double *A, szint lda, double *B, szint ldb)
  /* Solve AX=B, A is n-by-n, B is n-by-nrhs */
{
  szint info;
  szint *ipiv=new szint[n];
  DGESV(&n,&nrhs,A,&lda,ipiv,B,&ldb,&info);
  delete[] ipiv;
  return info;
}

szint cdgetrf(szint m, szint n, double *A, szint lda, szint *ipiv)
  /* Compute LU factorization of matrix A, m-by-n
     ipiv pivot indices, of size min(m,n) */
{
  szint info;
  DGETRF(&m,&n,A,&lda,ipiv,&info);
  return info;
}

szint cdgetrs(char trans, szint n, szint nrhs, double *A, szint lda,
            szint *ipiv, double *B, szint ldb)
  /* Solve AX=B for matrix X using LU from cdgetrf
     A is n-by-n, B is n-by-nrhs
     trans = 'N' (no transpose) or 'T' (transpose)
     Note: A is the output from cdgetrf, not the original matrix */
{
  szint info;
  DGETRS(&trans,&n,&nrhs,A,&lda,ipiv,B,&ldb,&info);
  return info;
}

szint cdgetri(szint n, double *A, szint lda, szint *ipiv)
  /* Compute inverse of matrix A using LU from cdgetrf
     Note: A is the output from cdgetrf, not the original matrix
     On exit, A is the inverse of the original matrix */
{
  szint info;
  szint lwork=16*n;
  double *work=new double[lwork];
  DGETRI(&n,A,&lda,ipiv,work,&lwork,&info);
  delete[] work;
  return info;
}

szint cdpotrf(char uplo, szint n, double *A, szint lda)
  /* Compute Cholesky factorization of SPD matrix A, n-by-n
     uplo = 'U' or 'L' for upper or lower triangular */
{
  szint info;
  DPOTRF(&uplo,&n,A,&lda,&info);
  return info;
}

szint cdpotrs(char uplo, szint n, szint nrhs, double *A, szint lda, double *B, szint ldb)
  /* Solve AX=B for SPD matrix A using Cholesky from cdpotrf
     A is n-by-n, B is n-by-nrhs, uplo as in cdpotfs
     Note: A is the output from cdpotfs, not the original matrix */
{
  szint info;
  DPOTRS(&uplo,&n,&nrhs,A,&lda,B,&ldb,&info);
  return info;
}

szint cdpotri(char uplo, szint n, double *A, szint lda)
  /* Compute inverse of matrix A using Cholesky from cdpotrf
     Note: A is the output from cdpotrf, not the original matrix
     On exit, A is the inverse of the original matrix */
{
  szint info;
  DPOTRI(&uplo,&n,A,&lda,&info);
  return info;
}

szint cdsyev(char jobz, char uplo, szint n, double *A, szint lda, double *w)
  /* Compute all eigenvalues and, optionally, eigenvectors of n-by-n real
     ymmetric matrix A. jobz = 'N' (do not compute eigenvectors) or 'V'
     (compute eigenvectors). Eigenvalues in w, eigenvectors in A
   */
{
  szint info;
  szint lwork=3*n;
  double *work=new double[lwork];
  DSYEV(&jobz,&uplo,&n,A,&lda,w,work,&lwork,&info);
  delete[] work;
  return info;
}

// szint cdgeev(char jobvl, char jobvr, szint n, double *A, szint lda, double *wr,
//            double *wi, double *vl, szint ldvl, double *vr, szint ldvr)
//   /* Compute eigenvalues and, optionally, left and/or right eigenvectors,
//      of n-by-n real nonsymmetric matrix A. jobvl = 'N' (do not compute
//      left eigenvectors) or 'V' (compute eigenvectors), similar for jobvr.
//      Eigenvalues in wr (real part) and wi (imaginary part). Overwrites A.
//    */
// {
//   szint info;
//   szint lwork=4*n;
//   double *work=new double[lwork];
//   DGEEV(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
//   delete[] work;
//   return info;
// }

void cdgemm(char transa, char transb, szint m, szint n, szint k, double alpha,
            double *A, szint lda, double *B, szint ldb, double beta, double *C, szint ldc)
  /* Matrix-matrix multiplication C = alpha op(A) op(B) + beta C
     A is m-by-k, B is k-by-n, C is m-by-n
     op depends on transa/transb = 'N' (no transpose) or 'T' (transpose) */
{
  DGEMM(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

void cdgemv(char trans, szint m, szint n, double alpha, double *A, szint lda,
            double *X, szint incx, double beta, double *Y, szint incy)
  /* Matrix-vector multiplication y = alpha op(A) x + beta y
     A is m-by-n, trans as in cdgemm */
{
  DGEMV(&trans,&m,&n,&alpha,A,&lda,X,&incx,&beta,Y,&incy);
}

void cdger(szint m, szint n, double alpha, double *x, szint incx, 
           double *y, szint incy, double *A, szint lda)
  /* Outer product A = alpha x y^T + A
     m elements in x, n elements in y, A is m-by-n */
{
  DGER(&m,&n,&alpha,x,&incx,y,&incy,A,&lda);
}

void cdaxpy(szint m, double alpha, double *x, szint incx, double *y, szint incy)
  /* y = y + alpha*x, with m elements in x,y */
{
  DAXPY(&m,&alpha,x,&incx,y,&incy);
}

void cdcopy(szint m, double *dx, szint incx, double *dy, szint incy)
  /* dy = dx, with m elements in dx,dy */
{
  DCOPY(&m,dx,&incx,dy,&incy);
}

void cdscal(szint m, double alpha, double *x, szint incx)
  /* x = alpha*x, with m elements in x */
{
  DSCAL(&m,&alpha,x,&incx);
}

double cdnrm2(szint n, double *x, szint incx)
  /* returns 2-norm of vector x with n elements */
{
  return DNRM2(&n,x,&incx);
}

//////////////////// SINGLE PRECISION /////////////////////////////////

void cslaswp(szint n, float *A, szint lda, szint k1, szint k2,
             szint *ipiv, szint incx)
  /* Perform row interchanges on A, from k1 to k2 in ipiv */
{
  SLASWP(&n,A,&lda,&k1,&k2,ipiv,&incx);
}

void cstrsm(char side, char uplo, char transa, char diag, szint m, szint n,
            float alpha, float *A, szint lda, float *B, szint ldb)
  /* Solve AX=alpha*B or XA=alpha*B for X with A triangular */
{
  STRSM(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}

void cstpsv(char uplo, char trans, char diag, szint n, float *Ap, float *x, szint incx)
  /* Solve Ax=b, A is n-by-n triangular */ 
{
  STPSV(&uplo,&trans,&diag,&n,Ap,x,&incx);
}

szint csgesv(szint n, szint nrhs, float *A, szint lda, float *B, szint ldb)
  /* Solve AX=B, A is n-by-n, B is n-by-nrhs */
{
  szint info;
  szint *ipiv=new szint[n];
  SGESV(&n,&nrhs,A,&lda,ipiv,B,&ldb,&info);
  delete[] ipiv;
  return info;
}

szint csgetrf(szint m, szint n, float *A, szint lda, szint *ipiv)
  /* Compute LU factorization of matrix A, m-by-n
     ipiv pivot indices, of size min(m,n) */
{
  szint info;
  SGETRF(&m,&n,A,&lda,ipiv,&info);
  return info;
}

szint csgetrs(char trans, szint n, szint nrhs, float *A, szint lda,
            szint *ipiv, float *B, szint ldb)
  /* Solve AX=B for matrix X using LU from csgetrf
     A is n-by-n, B is n-by-nrhs
     trans = 'N' (no transpose) or 'T' (transpose)
     Note: A is the output from csgetrf, not the original matrix */
{
  szint info;
  SGETRS(&trans,&n,&nrhs,A,&lda,ipiv,B,&ldb,&info);
  return info;
}

szint csgetri(szint n, float *A, szint lda, szint *ipiv)
  /* Compute inverse of matrix A using LU from csgetrf
     Note: A is the output from csgetrf, not the original matrix
     On exit, A is the inverse of the original matrix */
{
  szint info;
  szint lwork=16*n;
  float *work=new float[lwork];
  SGETRI(&n,A,&lda,ipiv,work,&lwork,&info);
  delete[] work;
  return info;
}

szint cspotrf(char uplo, szint n, float *A, szint lda)
  /* Compute Cholesky factorization of SPD matrix A, n-by-n
     uplo = 'U' or 'L' for upper or lower triangular */
{
  szint info;
  SPOTRF(&uplo,&n,A,&lda,&info);
  return info;
}

szint cspotrs(char uplo, szint n, szint nrhs, float *A, szint lda, float *B, szint ldb)
  /* Solve AX=B for SPD matrix A using Cholesky from cspotrf
     A is n-by-n, B is n-by-nrhs, uplo as in cspotfs
     Note: A is the output from cspotfs, not the original matrix */
{
  szint info;
  SPOTRS(&uplo,&n,&nrhs,A,&lda,B,&ldb,&info);
  return info;
}

szint cspotri(char uplo, szint n, float *A, szint lda)
  /* Compute inverse of matrix A using Cholesky from cspotrf
     Note: A is the output from cspotrf, not the original matrix
     On exit, A is the inverse of the original matrix */
{
  szint info;
  SPOTRI(&uplo,&n,A,&lda,&info);
  return info;
}

szint cssyev(char jobz, char uplo, szint n, float *A, szint lda, float *w)
  /* Compute all eigenvalues and, optionally, eigenvectors of n-by-n real
     ymmetric matrix A. jobz = 'N' (do not compute eigenvectors) or 'V'
     (compute eigenvectors). Eigenvalues in w, eigenvectors in A
   */
{
  szint info;
  szint lwork=3*n;
  float *work=new float[lwork];
  SSYEV(&jobz,&uplo,&n,A,&lda,w,work,&lwork,&info);
  delete[] work;
  return info;
}

// szint csgeev(char jobvl, char jobvr, szint n, float *A, szint lda, float *wr,
//            float *wi, float *vl, szint ldvl, float *vr, szint ldvr)
//   /* Compute eigenvalues and, optionally, left and/or right eigenvectors,
//      of n-by-n real nonsymmetric matrix A. jobvl = 'N' (do not compute
//      left eigenvectors) or 'V' (compute eigenvectors), similar for jobvr.
//      Eigenvalues in wr (real part) and wi (imaginary part). Overwrites A.
//    */
// {
//   szint info;
//   szint lwork=4*n;
//   float *work=new float[lwork];
//   SGEEV(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
//   delete[] work;
//   return info;
// }

void csgemm(char transa, char transb, szint m, szint n, szint k, float alpha,
            float *A, szint lda, float *B, szint ldb, float beta, float *C, szint ldc)
  /* Matrix-matrix multiplication C = alpha op(A) op(B) + beta C
     A is m-by-k, B is k-by-n, C is m-by-n
     op depends on transa/transb = 'N' (no transpose) or 'T' (transpose) */
{
  SGEMM(&transa,&transb,&m,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

void csgemv(char trans, szint m, szint n, float alpha, float *A, szint lda,
            float *X, szint incx, float beta, float *Y, szint incy)
  /* Matrix-vector multiplication y = alpha op(A) x + beta y
     A is m-by-n, trans as in csgemm */
{
  SGEMV(&trans,&m,&n,&alpha,A,&lda,X,&incx,&beta,Y,&incy);
}

void csger(szint m, szint n, float alpha, float *x, szint incx, 
           float *y, szint incy, float *A, szint lda)
  /* Outer product A = alpha x y^T + A
     m elements in x, n elements in y, A is m-by-n */
{
  SGER(&m,&n,&alpha,x,&incx,y,&incy,A,&lda);
}

void csaxpy(szint m, float alpha, float *x, szint incx, float *y, szint incy)
  /* y = y + alpha*x, with m elements in x,y */
{
  SAXPY(&m,&alpha,x,&incx,y,&incy);
}

void cscopy(szint m, float *dx, szint incx, float *dy, szint incy)
  /* dy = dx, with m elements in dx,dy */
{
  SCOPY(&m,dx,&incx,dy,&incy);
}

void csscal(szint m, float alpha, float *x, szint incx)
  /* x = alpha*x, with m elements in x */
{
  SSCAL(&m,&alpha,x,&incx);
}

float csnrm2(szint n, float *x, szint incx)
  /* returns 2-norm of vector x with n elements */
{
  return SNRM2(&n,x,&incx);
}
