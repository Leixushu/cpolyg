// FORTWRAP.H
//
// C-wrappers for BLAS and LAPACK fortran routines
//
// See doc for quick reference guides
// See fortwrap.doc for Fortran comments
// See C-comments below for short descriptions

#ifndef __FORTWRAP
#define __FORTWRAP

#include "fortblas.h"

/*
  General variable names:
 
  - lda, ldb, etc are leading dimensions of matrices A, B, etc
  - incx, incy, etc are increments of elements in x, y, etc
*/

//////////////////// DOUBLE PRECISION /////////////////////////////////

void cdlaswp(szint n, double *A, szint lda, szint k1, szint k2,
             szint *ipiv, szint incx);
  /* Perform row interchanges on A, from k1 to k2 in ipiv */

void cdtrsm(char side, char uplo, char transa, char diag, szint m, szint n,
            double alpha, double *A, szint lda, double *B, szint ldb);
  /* Solve AX=alpha*B or XA=alpha*B for X with A triangular */

void cdtpsv(char uplo, char trans, char diag, szint n, double *Ap, double *x, szint incx);
  /* Solve Ax=b, A is n-by-n triangular */ 

szint cdgesv(szint n, szint nrhs, double *A, szint lda, double *B, szint ldb);
  /* Solve AX=B, A is n-by-n, B is n-by-nrhs */

szint cdgetrf(szint m, szint n, double *A, szint lda, szint *ipiv);
  /* Compute LU factorization of matrix A, m-by-n
     ipiv pivot indices, of size min(m,n) */

szint cdgetrs(char trans, szint n, szint nrhs, double *A, szint lda,
            szint *ipiv, double *B, szint ldb);
  /* Solve AX=B for matrix A using LU from cdgetrf
     A is n-by-n, B is n-by-nrhs
     trans = 'N' (no transpose) or 'T' (transpose)
     Note: A is the output from cdgetrf, not the original matrix */

szint cdgetri(szint n, double *A, szint lda, szint *ipiv);
  /* Compute inverse of matrix A using LU from cdgetrf
     Note: A is the output from cdgetrf, not the original matrix
     On exit, A is the inverse of the original matrix */

szint cdpotrf(char uplo, szint n, double *A, szint lda);
  /* Compute Cholesky factorization of SPD matrix A, n-by-n
     uplo = 'U' or 'L' for upper or lower triangular */

szint cdpotrs(char uplo, szint n, szint nrhs, double *A, szint lda, double *B, szint ldb);
  /* Solve AX=B for SPD matrix A using Cholesky from cdpotrf
     A is n-by-n, B is n-by-nrhs, uplo as in cdpotfs
     Note: A is the output from cdpotfs, not the original matrix */

szint cdpotri(char uplo, szint n, double *A, szint lda);
  /* Compute inverse of matrix A using Cholesky from cdpotrf
     Note: A is the output from cdpotrf, not the original matrix
     On exit, A is the inverse of the original matrix */

szint cdsyev(char jobz, char uplo, szint n, double *A, szint lda, double *w);
  /* Compute all eigenvalues and, optionally, eigenvectors of n-by-n real
     ymmetric matrix A. jobz = 'N' (do not compute eigenvectors) or 'V'
     (compute eigenvectors). Eigenvalues in w, eigenvectors in A
   */

void cdgemm(char transa, char transb, szint m, szint n, szint k, double alpha,
            double *A, szint lda, double *B, szint ldb, double beta, double *C, szint ldc);
  /* Matrix-matrix multiplication C = alpha op(A) op(B) + beta C
     A is m-by-k, B is k-by-n, C is m-by-n
     op depends on transa/transb = 'N' (no transpose) or 'T' (transpose) */

void cdgemv(char trans, szint m, szint n, double alpha, double *A, szint lda,
            double *X, szint incx, double beta, double *Y, szint incy);
  /* Matrix-vector multiplication y = alpha op(A) x + beta y
     A is m-by-n, trans as in cdgemm */

void cdger(szint m, szint n, double alpha, double *x, szint incx, 
           double *y, szint incy, double *A, szint lda);
  /* Outer product A = alpha x y^T + A
     m elements in x, n elements in y, A is m-by-n */

void cdaxpy(szint m, double alpha, double *x, szint incx, double *y, szint incy);
  /* y = y + alpha*x, with m elements in x,y */

void cdcopy(szint m, double *dx, szint incx, double *dy, szint incy);
  /* dy = dx, with m elements in dx,dy */

void cdscal(szint m, double alpha, double *x, szint incx);
  /* x = alpha*x, with m elements in x */

double cdnrm2(szint n, double *x, szint incx);
  /* returns 2-norm of vector x with n elements */

//////////////////// SINGLE PRECISION /////////////////////////////////

void cslaswp(szint n, float *A, szint lda, szint k1, szint k2,
             szint *ipiv, szint incx);
  /* Perform row interchanges on A, from k1 to k2 in ipiv */

void cstrsm(char side, char uplo, char transa, char diag, szint m, szint n,
            float alpha, float *A, szint lda, float *B, szint ldb);
  /* Solve AX=alpha*B or XA=alpha*B for X with A triangular */

void cstpsv(char uplo, char trans, char diag, szint n, float *Ap, float *x, szint incx);
  /* Solve Ax=b, A is n-by-n triangular */ 

szint csgesv(szint n, szint nrhs, float *A, szint lda, float *B, szint ldb);
  /* Solve AX=B, A is n-by-n, B is n-by-nrhs */

szint csgetrf(szint m, szint n, float *A, szint lda, szint *ipiv);
  /* Compute LU factorization of matrix A, m-by-n
     ipiv pivot indices, of size min(m,n) */

szint csgetrs(char trans, szint n, szint nrhs, float *A, szint lda,
            szint *ipiv, float *B, szint ldb);
  /* Solve AX=B for matrix A using LU from csgetrf
     A is n-by-n, B is n-by-nrhs
     trans = 'N' (no transpose) or 'T' (transpose)
     Note: A is the output from csgetrf, not the original matrix */

szint csgetri(szint n, float *A, szint lda, szint *ipiv);
  /* Compute inverse of matrix A using LU from csgetrf
     Note: A is the output from csgetrf, not the original matrix
     On exit, A is the inverse of the original matrix */

szint cspotrf(char uplo, szint n, float *A, szint lda);
  /* Compute Cholesky factorization of SPD matrix A, n-by-n
     uplo = 'U' or 'L' for upper or lower triangular */

szint cspotrs(char uplo, szint n, szint nrhs, float *A, szint lda, float *B, szint ldb);
  /* Solve AX=B for SPD matrix A using Cholesky from cspotrf
     A is n-by-n, B is n-by-nrhs, uplo as in cspotfs
     Note: A is the output from cspotfs, not the original matrix */

szint cspotri(char uplo, szint n, float *A, szint lda);
  /* Compute inverse of matrix A using Cholesky from cspotrf
     Note: A is the output from cspotrf, not the original matrix
     On exit, A is the inverse of the original matrix */

szint cssyev(char jobz, char uplo, szint n, float *A, szint lda, float *w);
  /* Compute all eigenvalues and, optionally, eigenvectors of n-by-n real
     ymmetric matrix A. jobz = 'N' (do not compute eigenvectors) or 'V'
     (compute eigenvectors). Eigenvalues in w, eigenvectors in A
   */

void csgemm(char transa, char transb, szint m, szint n, szint k, float alpha,
            float *A, szint lda, float *B, szint ldb, float beta, float *C, szint ldc);
  /* Matrix-matrix multiplication C = alpha op(A) op(B) + beta C
     A is m-by-k, B is k-by-n, C is m-by-n
     op depends on transa/transb = 'N' (no transpose) or 'T' (transpose) */

void csgemv(char trans, szint m, szint n, float alpha, float *A, szint lda,
            float *X, szint incx, float beta, float *Y, szint incy);
  /* Matrix-vector multiplication y = alpha op(A) x + beta y
     A is m-by-n, trans as in csgemm */

void csger(szint m, szint n, float alpha, float *x, szint incx, 
           float *y, szint incy, float *A, szint lda);
  /* Outer product A = alpha x y^T + A
     m elements in x, n elements in y, A is m-by-n */

void csaxpy(szint m, float alpha, float *x, szint incx, float *y, szint incy);
  /* y = y + alpha*x, with m elements in x,y */

void cscopy(szint m, float *dx, szint incx, float *dy, szint incy);
  /* dy = dx, with m elements in dx,dy */

void csscal(szint m, float alpha, float *x, szint incx);
  /* x = alpha*x, with m elements in x */

float csnrm2(szint n, float *x, szint incx);
  /* returns 2-norm of vector x with n elements */

#endif
