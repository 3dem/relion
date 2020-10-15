#ifndef __LIBGRAVIS_LAPACK_H__
#define __LIBGRAVIS_LAPACK_H__
/******************************************************************************
 **        Title: lapack.h
 **  Description: Connect our matrix classes to lapack.
 **
 **       Author: Brian Amberg, 2007
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#include "Exception.h"
#include "t2Matrix.h"
#include "t3Matrix.h"
#include "t4Matrix.h"

// Lapack Header
extern "C" {
  void sgesvd_(const char& jobu, const char& jobvt, const int& m, const int& n, float* a, const int& lda, float* s, float* u, const int& ldu, float* vt, const int& ldvt, float* work, const int& lwork, int& info );
  void dgesvd_(const char& jobu, const char& jobvt, const int& m, const int& n, double* a, const int& lda, double* s, double* u, const int& ldu, double* vt, const int& ldvt, double* work, const int& lwork, int& info );
}


namespace gravis
{
  /**
   *Use lapack to calculate an svd of a 2x2 matrix
   **/
  void svd(f2Matrix& U, f2Vector& S, f2Matrix& VT, const f2Matrix& A)
  {
    f2Matrix _A(A);
    float WORK[16];
    int INFO;
    sgesvd_('A', 'A', 2, 2, &(_A[0]), 2, &(S[0]), &(U[0]), 2, &(VT[0]), 2,  WORK, 16, INFO);
    if (INFO < 0) GRAVIS_THROW2(gravis::Exception, "The i'th argument had an invalid value.");
    if (INFO > 0) GRAVIS_THROW2(gravis::Exception, "SBDSQR did not converge to zero.");
  }
  /**
   *Use lapack to calculate an svd of a 2x2 matrix
   **/
  void svd(d2Matrix& U, d2Vector& S, d2Matrix& VT, const d2Matrix& A)
  {
    d2Matrix _A(A);
    double WORK[16];
    int INFO;
    dgesvd_('A', 'A', 2, 2, &(_A[0]), 2, &(S[0]), &(U[0]), 2, &(VT[0]), 2,  WORK, 16, INFO);
    if (INFO < 0) GRAVIS_THROW2(gravis::Exception, "The i'th argument had an invalid value.");
    if (INFO > 0) GRAVIS_THROW2(gravis::Exception, "SBDSQR did not converge to zero.");
  }

  /**
   *Use lapack to calculate an svd of a 3x3 matrix
   **/
  void svd(f3Matrix& U, f3Vector& S, f3Matrix& VT, const f3Matrix& A)
  {
    f3Matrix _A(A);
    float WORK[32];
    int INFO;
    sgesvd_('A', 'A', 3, 3, &(_A[0]), 3, &(S[0]), &(U[0]), 3, &(VT[0]), 3,  WORK, 32, INFO);
    if (INFO < 0) GRAVIS_THROW2(gravis::Exception, "The i'th argument had an invalid value.");
    if (INFO > 0) GRAVIS_THROW2(gravis::Exception, "SBDSQR did not converge to zero.");
  }
  /**
   *Use lapack to calculate an svd of a 3x3 matrix
   **/
  void svd(d3Matrix& U, d3Vector& S, d3Matrix& VT, const d3Matrix& A)
  {
    d3Matrix _A(A);
    double WORK[32];
    int INFO;
    dgesvd_('A', 'A', 3, 3, &(_A[0]), 3, &(S[0]), &(U[0]), 3, &(VT[0]), 3,  WORK, 32, INFO);
    if (INFO < 0) GRAVIS_THROW2(gravis::Exception, "The i'th argument had an invalid value.");
    if (INFO > 0) GRAVIS_THROW2(gravis::Exception, "SBDSQR did not converge to zero.");
  }

  /**
   *Use lapack to calculate an svd of a 4x4 matrix
   **/
  void svd(f4Matrix& U, f4Vector& S, f4Matrix& VT, const f4Matrix& A)
  {
    f4Matrix _A(A);
    float WORK[64];
    int INFO;
    sgesvd_('A', 'A', 4, 4, &(_A[0]), 4, &(S[0]), &(U[0]), 4, &(VT[0]), 4,  WORK, 64, INFO);
    if (INFO < 0) GRAVIS_THROW2(gravis::Exception, "The i'th argument had an invalid value.");
    if (INFO > 0) GRAVIS_THROW2(gravis::Exception, "SBDSQR did not converge to zero.");
  }
  /**
   *Use lapack to calculate an svd of a 4x4 matrix
   **/
  void svd(d4Matrix& U, d4Vector& S, d4Matrix& VT, const d4Matrix& A)
  {
    d4Matrix _A(A);
    double WORK[64];
    int INFO;
    dgesvd_('A', 'A', 4, 4, &(_A[0]), 4, &(S[0]), &(U[0]), 4, &(VT[0]), 4,  WORK, 64, INFO);
    if (INFO < 0) GRAVIS_THROW2(gravis::Exception, "The i'th argument had an invalid value.");
    if (INFO > 0) GRAVIS_THROW2(gravis::Exception, "SBDSQR did not converge to zero.");
  }

  int rank(const f2Matrix& A, const float accuracy = 1e-10)
  {
    f2Matrix U;
    f2Vector S;
    f2Matrix VT;
    svd(U, S, VT, A);
    int r = 0;
    while (r<2 && (S[r] >= accuracy || S[r] <= -accuracy)) ++r;
    return r;
  }

  int rank(const d2Matrix& A, const double accuracy = 1e-10)
  {
    d2Matrix U;
    d2Vector S;
    d2Matrix VT;
    svd(U, S, VT, A);
    int r = 0;
    while (r<2 && (S[r] >= accuracy || S[r] <= -accuracy)) ++r;
    return r;
  }

  int rank(const f3Matrix& A, const float accuracy = 1e-10)
  {
    f3Matrix U;
    f3Vector S;
    f3Matrix VT;
    svd(U, S, VT, A);
    int r = 0;
    while (r<3 && (S[r] >= accuracy || S[r] <= -accuracy)) ++r;
    return r;
  }

  int rank(const d3Matrix& A, const double accuracy = 1e-10)
  {
    d3Matrix U;
    d3Vector S;
    d3Matrix VT;
    svd(U, S, VT, A);
    int r = 0;
    while (r<3 && (S[r] >= accuracy || S[r] <= -accuracy)) ++r;
    return r;
  }

  int rank(const f4Matrix& A, const float accuracy = 1e-10)
  {
    f4Matrix U;
    f4Vector S;
    f4Matrix VT;
    svd(U, S, VT, A);
    int r = 0;
    while (r<4 && (S[r] >= accuracy || S[r] <= -accuracy)) ++r;
    return r;
  }

  int rank(const d4Matrix& A, const double accuracy = 1e-10)
  {
    d4Matrix U;
    d4Vector S;
    d4Matrix VT;
    svd(U, S, VT, A);
    int r = 0;
    while (r<4 && (S[r] >= accuracy || S[r] <= -accuracy)) ++r;
    return r;
  }
}
#endif
