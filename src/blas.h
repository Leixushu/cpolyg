#ifndef __BLAS
#define __BLAS

#include <cstdlib>
#include <cstddef>

#ifdef MATLAB_MEX_FILE
#define szint ptrdiff_t
#else
#define szint int
#endif

#if defined(_WIN32)
#define ADD_UNDERSCORE(x) x
#else
#define ADD_UNDERSCORE(x) x ## _
#endif

#define XERBLA ADD_UNDERSCORE(xerbla)

#define DLASWP ADD_UNDERSCORE(dlaswp)
#define DTRSM  ADD_UNDERSCORE(dtrsm)
#define DTPSV  ADD_UNDERSCORE(dtpsv)
#define DGESV  ADD_UNDERSCORE(dgesv)
#define DGETRF ADD_UNDERSCORE(dgetrf)
#define DGETRS ADD_UNDERSCORE(dgetrs)
#define DGETRI ADD_UNDERSCORE(dgetri)
#define DPOTRF ADD_UNDERSCORE(dpotrf)
#define DPOTRS ADD_UNDERSCORE(dpotrs)
#define DPOTRI ADD_UNDERSCORE(dpotri)
#define DSYEV  ADD_UNDERSCORE(dsyev)
#define DGEEV  ADD_UNDERSCORE(dgeev)
#define DGEMM  ADD_UNDERSCORE(dgemm)
#define DGEMV  ADD_UNDERSCORE(dgemv)
#define DGER   ADD_UNDERSCORE(dger)
#define DAXPY  ADD_UNDERSCORE(daxpy)
#define DCOPY  ADD_UNDERSCORE(dcopy)
#define DSCAL  ADD_UNDERSCORE(dscal)
#define DNRM2  ADD_UNDERSCORE(dnrm2)

#define SLASWP ADD_UNDERSCORE(slaswp)
#define STRSM  ADD_UNDERSCORE(strsm)
#define STPSV  ADD_UNDERSCORE(stpsv)
#define SGESV  ADD_UNDERSCORE(sgesv)
#define SGETRF ADD_UNDERSCORE(sgetrf)
#define SGETRS ADD_UNDERSCORE(sgetrs)
#define SGETRI ADD_UNDERSCORE(sgetri)
#define SPOTRF ADD_UNDERSCORE(spotrf)
#define SPOTRS ADD_UNDERSCORE(spotrs)
#define SPOTRI ADD_UNDERSCORE(spotri)
#define SSYEV  ADD_UNDERSCORE(ssyev)
#define SGEEV  ADD_UNDERSCORE(sgeev)
#define SGEMM  ADD_UNDERSCORE(sgemm)
#define SGEMV  ADD_UNDERSCORE(sgemv)
#define SGER   ADD_UNDERSCORE(sger)
#define SAXPY  ADD_UNDERSCORE(saxpy)
#define SCOPY  ADD_UNDERSCORE(scopy)
#define SSCAL  ADD_UNDERSCORE(sscal)
#define SNRM2  ADD_UNDERSCORE(snrm2)

/* extern "C" void XERBLA()  { } */
extern "C" {
  void DLASWP(szint*,double*,szint*,szint*,szint*,szint*,szint*);
  void DTRSM(char*,char*,char*,char*,szint*,szint*,double*,double*,szint*,
             double*,szint*);
  void DTPSV(char*,char*,char*,szint*,double*,double*,szint*);
  void DGESV(szint*,szint*,double*,szint*,szint*,double*,szint*,szint*);
  void DGETRF(szint*,szint*,double*,szint*,szint*,szint*);
  void DGETRS(char*,szint*,szint*,double*,szint*,szint*,double*,szint*,szint*);
  void DGETRI(szint*,double*,szint*,szint*,double*,szint*,szint*);
  void DPOTRF(char*,szint*,double*,szint*,szint*);
  void DPOTRS(char*,szint*,szint*,double*,szint*,double*,szint*,szint*);
  void DPOTRI(char*,szint*,double*,szint*,szint*);
  void DSYEV(char*,char*,szint*,double*,szint*,double*,double*,szint*,szint*);
  void DGEEV(char*,char*,szint*,double*,szint*,double*,double*,
             double*,szint*,double*,szint*,double*,szint*,szint*);
  void DGEMM(char*,char*,szint*,szint*,szint*,double*,double*,szint*,
             double*,szint*,double*,double*,szint*);
  void DGEMV(char*,szint*,szint*,double*,double*,szint*,double*,szint*,double*,double*,szint*);
  void DGER(szint*,szint*,double*,double*,szint*,double*,szint*,double*,szint*);
  void DAXPY(szint*,double*,double*,szint*,double*,szint*);
  void DCOPY(szint*,double*,szint*,double*,szint*);
  void DSCAL(szint*,double*,double*,szint*);
  double DNRM2(szint*,double*,szint*);

  void SLASWP(szint*,float*,szint*,szint*,szint*,szint*,szint*);
  void STRSM(char*,char*,char*,char*,szint*,szint*,float*,float*,szint*,
             float*,szint*);
  void STPSV(char*,char*,char*,szint*,float*,float*,szint*);
  void SGESV(szint*,szint*,float*,szint*,szint*,float*,szint*,szint*);
  void SGETRF(szint*,szint*,float*,szint*,szint*,szint*);
  void SGETRS(char*,szint*,szint*,float*,szint*,szint*,float*,szint*,szint*);
  void SGETRI(szint*,float*,szint*,szint*,float*,szint*,szint*);
  void SPOTRF(char*,szint*,float*,szint*,szint*);
  void SPOTRS(char*,szint*,szint*,float*,szint*,float*,szint*,szint*);
  void SPOTRI(char*,szint*,float*,szint*,szint*);
  void SSYEV(char*,char*,szint*,float*,szint*,float*,float*,szint*,szint*);
  void SGEEV(char*,char*,szint*,float*,szint*,float*,float*,
             float*,szint*,float*,szint*,float*,szint*,szint*);
  void SGEMM(char*,char*,szint*,szint*,szint*,float*,float*,szint*,
             float*,szint*,float*,float*,szint*);
  void SGEMV(char*,szint*,szint*,float*,float*,szint*,float*,szint*,float*,float*,szint*);
  void SGER(szint*,szint*,float*,float*,szint*,float*,szint*,float*,szint*);
  void SAXPY(szint*,float*,float*,szint*,float*,szint*);
  void SCOPY(szint*,float*,szint*,float*,szint*);
  void SSCAL(szint*,float*,float*,szint*);
  float SNRM2(szint*,float*,szint*);
}

#endif
