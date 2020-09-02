/**
 * Included multiple times from matrix_blas.hpp for different combinations of float, double varmatrix and matrixview
 *
 * Never include directly
 **/


#ifdef __GRAVIS__MATRIX__BLAS__DATATYPE__SINGLE__

#define __GMBD_REAL  float
#define __GMBD_xGEMV sgemv_
#define __GMBD_xNRM2 snrm2_
#define __GMBD_xSCAL sscal_
#define __GMBD_xAXPY saxpy_

#define __GMBD_xGESVD sgesvd_
#define __GMBD_xGESDD sgesdd_
#define __GMBD_xDOT   sdot_

#else
#ifdef __GRAVIS__MATRIX__BLAS__DATATYPE__DOUBLE__

#define __GMBD_REAL  double
#define __GMBD_xGEMV dgemv_
#define __GMBD_xNRM2 dnrm2_
#define __GMBD_xSCAL dscal_
#define __GMBD_xAXPY daxpy_

#define __GMBD_xGESVD dgesvd_
#define __GMBD_xGESDD dgesdd_
#define __GMBD_xDOT ddot_

#else
#error( "Never include directly, this is included only from within matrix_blas.hpp" )
#endif
#endif

namespace reference
{
  //#include "tVarMatrix_blas_reference.h"
}

// Blas Header
extern "C" {
  void        __GMBD_xGEMV(const char* const trans, const size_t& m, const size_t& n,
                           const __GMBD_REAL& alpha, const __GMBD_REAL* const M, const size_t& m1, const __GMBD_REAL* const x, const size_t& xs,
                           const __GMBD_REAL& beta, __GMBD_REAL* const v, const size_t& vs);
  __GMBD_REAL __GMBD_xNRM2(const size_t& n, const __GMBD_REAL* const x, const size_t& inc);
  void        __GMBD_xSCAL(const size_t& n, const __GMBD_REAL& alpha, __GMBD_REAL* const x, const size_t& inc);
  void        __GMBD_xAXPY(const size_t& n, const __GMBD_REAL& alpha, const __GMBD_REAL* const x, const size_t& incx, __GMBD_REAL* const y, const size_t& incy);
  __GMBD_REAL __GMBD_xDOT(const size_t& n, const __GMBD_REAL* dx, const size_t& incx, const __GMBD_REAL* dy, const size_t& incy);
}

/// Lapack Header
extern "C" {
  void __GMBD_xGESVD(const char& jobu, const char& jobvt, const int& m, const int& n, __GMBD_REAL* a, const int& lda, __GMBD_REAL* s, __GMBD_REAL* u, const int& ldu, __GMBD_REAL* vt, const int& ldvt, __GMBD_REAL* work, const int& lwork, int& info );
  void __GMBD_xGESDD(const char& jobz, const int& m, const int& n, __GMBD_REAL* a, const int& lda, __GMBD_REAL* s, __GMBD_REAL* u, const int& ldu, __GMBD_REAL* vt, const int& ldvt, __GMBD_REAL* work, const int& lwork, int* iwork, int& info );
}

/**
 * Inplace SVD for small matrices.
 * Replaces the input matrix A with its left eigenvectors U
 **/
inline static
void svd_inplace_u(tMatrixView<__GMBD_REAL> &IN_A_OUT_U, tVectorView<__GMBD_REAL> &S, tMatrixView<__GMBD_REAL> &VT)
{
  int info;
  tVarVector<__GMBD_REAL> work(1);
  __GMBD_xGESVD('O', 'A', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, -1, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESVD('O', 'A', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, work.h, info);
}

/**
 * Inplace SVD for small matrices
 * Replaces the input matrix A with its left eigenvectors U
 **/
inline static
void svd_inplace_u(tVarMatrix<__GMBD_REAL> &IN_A_OUT_U, tVarVector<__GMBD_REAL> &S, tVarMatrix<__GMBD_REAL> &VT)
{
  int info;
  tVarVector<__GMBD_REAL> work(1);
  __GMBD_xGESVD('O', 'A', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, -1, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESVD('O', 'A', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, work.h, info);
}

/**
 * Inplace SVD for large matrices using a divide and conquer algorithm
 * Replaces the input matrix A with its left eigenvectors U
 **/
inline static
void svd_inplace_u_dc(tMatrixView<__GMBD_REAL> &IN_A_OUT_U, tVectorView<__GMBD_REAL> &S, tMatrixView<__GMBD_REAL> &VT)
{
  int info;
  tVarVector<__GMBD_REAL> work(1);
  tVarVector<int> iwork(8*std::min(IN_A_OUT_U.h, IN_A_OUT_U.w));
  __GMBD_xGESDD('O', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, -1, iwork.data, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESDD('O', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, work.h, iwork.data, info);
}

/**
 * Inplace SVD for large matrices using a divide and conquer algorithm
 * Replaces the input matrix A with its left eigenvectors U
 **/
inline static
void svd_inplace_u_dc(tVarMatrix<__GMBD_REAL> &IN_A_OUT_U, tVarVector<__GMBD_REAL> &S, tVarMatrix<__GMBD_REAL> &VT)
{
  int info;
  tVarVector<__GMBD_REAL> work(1);
  tVarVector<int> iwork(8*std::min(IN_A_OUT_U.h, IN_A_OUT_U.w));
  __GMBD_xGESDD('O', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, -1, iwork.data, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESDD('O', IN_A_OUT_U.h, IN_A_OUT_U.w, IN_A_OUT_U.data, IN_A_OUT_U.h, S.data, NULL, IN_A_OUT_U.h, VT.data, VT.h, work.data, work.h, iwork.data, info);
}


/**
 * SVD for small matrices
 **/
inline static
void svd(tMatrixView<__GMBD_REAL> &U, tVectorView<__GMBD_REAL> &S, tMatrixView<__GMBD_REAL> &VT, const tConstMatrixView<__GMBD_REAL> &A)
{
  int info;
  tVarMatrix<__GMBD_REAL> _A(A);
  tVarVector<__GMBD_REAL> work(1);
  __GMBD_xGESVD('A', 'A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, -1, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESVD('A', 'A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, work.h, info);
}

/**
 * SVD for small matrices
 **/
inline static
void svd(tVarMatrix<__GMBD_REAL> &U, tVarVector<__GMBD_REAL> &S, tVarMatrix<__GMBD_REAL> &VT, const tConstMatrixView<__GMBD_REAL> &A)
{
  int info;
  tVarMatrix<__GMBD_REAL> _A(A);
  tVarVector<__GMBD_REAL> work(1);
  __GMBD_xGESVD('A', 'A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, -1, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESVD('A', 'A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, work.h, info);
}

/**
 * SVD for large matrices using a divide and conquer algorithm
 **/
inline static
void svd_dc(tMatrixView<__GMBD_REAL> &U, tVectorView<__GMBD_REAL> &S, tMatrixView<__GMBD_REAL> &VT, const tConstMatrixView<__GMBD_REAL> &A)
{
  int info;
  tVarMatrix<__GMBD_REAL> _A(A);
  tVarVector<__GMBD_REAL> work(1);
  tVarVector<int> iwork(8*std::min(A.h, A.w));
  __GMBD_xGESDD('A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, -1, iwork.data, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESDD('A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, work.h, iwork.data, info);
}

/**
 * SVD for large matrices using a divide and conquer algorithm
 **/
inline static
void svd_dc(tVarMatrix<__GMBD_REAL> &U, tVarVector<__GMBD_REAL> &S, tVarMatrix<__GMBD_REAL> &VT, const tConstMatrixView<__GMBD_REAL> &A)
{
  int info;
  tVarMatrix<__GMBD_REAL> _A(A);
  tVarVector<__GMBD_REAL> work(1);
  tVarVector<int> iwork(8*std::min(A.h, A.w));
  __GMBD_xGESDD('A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, -1, iwork.data, info);
  if (info < 0) GRAVIS_THROW3(gravis::Exception, "The i'th argument had an invalid value.", StringFormat(info));
  if (info > 0) GRAVIS_THROW3(gravis::Exception, "SBDSQR did not converge to zero.",        StringFormat(info));
  work.resize(int(work[0]));
  __GMBD_xGESDD('A', A.h, A.w, _A.data, A.h, S.data, U.data, U.h, VT.data, VT.h, work.data, work.h, iwork.data, info);
}


//// Multiplications
// Not Transposed
/**
 * v = alpha*M*x + beta*v
 **/
inline static
void addmult(tVectorView<__GMBD_REAL> &v,
             const __GMBD_REAL& alpha, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x,
             const __GMBD_REAL& beta)
{
  GRAVIS_CHECK( v.size()  == M.h, "v and M are incompatible");
  GRAVIS_CHECK( x.size()  == M.w, "M and x are incompatible");
  if (M.h > 0)
    __GMBD_xGEMV("N", M.h, M.w, alpha, M.data, M.h, x.data, 1, beta, v.data, 1);
}
/**
 * v = alpha*M*x + beta*v
 * Will not resize v, as this would not make sense
 **/
inline static
void addmult(tVarVector<__GMBD_REAL> &v,
             const __GMBD_REAL& alpha, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x,
             const __GMBD_REAL& beta)
{
  tVectorView<__GMBD_REAL> vv(v);
  addmult(vv, alpha, M, x, beta);
}

/**
 * v = v+M*x
 **/
inline static
void addmult(tVectorView<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  addmult(v, __GMBD_REAL(1), M, x, __GMBD_REAL(1));
}
/**
 * v = v+M*x
 * Will not resize v, as this would not make sense
 **/
inline static
void addmult(tVarVector<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  addmult(v, __GMBD_REAL(1), M, x, __GMBD_REAL(1));
}

/**
 * v = a+M*x
 **/
inline static
void addmult(tVectorView<__GMBD_REAL> &v,
             const tConstVectorView<__GMBD_REAL> &a, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  v = a;
  addmult(v, M, x);
}

/**
 * v = a+M*x
 **/
inline static
void addmult(tVarVector<__GMBD_REAL> &v,
             const tConstVectorView<__GMBD_REAL> &a, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  v = a;
  addmult(v, M, x);
}
/**
 * v = alpha*M*x
 **/
inline static
void mult(tVectorView<__GMBD_REAL> &v,
          const __GMBD_REAL& alpha, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  ::gravis::matrix::clear(v);
  addmult(v, alpha, M, x, 1);
}
/**
 * v = alpha*M*x
 * Will not resize v, as this would not make sense
 **/
inline static
void mult(tVarVector<__GMBD_REAL> &v,
          const __GMBD_REAL& alpha, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  ::gravis::matrix::clear(v);
  addmult(v, alpha, M, x, 1);
}

/**
 * v = M*x
 **/
inline static
void mult(tVectorView<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  mult(v, 1, M, x);
}

/**
 * v = M*x
 **/
inline static
void mult(tVarVector<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  mult(v, 1, M, x);
}


// TRANSPOSED VERSIONS
/**
 * v = (alpha*x^T M)^T + beta*v
 **/
inline static
void addmult(tVectorView<__GMBD_REAL> &v,
             const __GMBD_REAL& alpha, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M,
             const __GMBD_REAL& beta)
{
  GRAVIS_CHECK( v.size() == M.w, "v and M are incompatible");
  GRAVIS_CHECK( x.size() == M.h, "M and x are incompatible");
  if (M.h > 0)
    __GMBD_xGEMV("T", M.h, M.w, alpha, M.data, M.h, x.data, 1, beta, v.data, 1);
}

/**
 * v = (alpha*x^T M)^T + beta*v
 **/
inline static
void addmult(tVarVector<__GMBD_REAL> &v,
             const __GMBD_REAL& alpha, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M,
             const __GMBD_REAL& beta)
{
  tVectorView<__GMBD_REAL> vv(v);
  addmult(vv, alpha, x, M, beta);
}
/**
 * v = v+M*x
 **/
inline static
void addmult(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  addmult(v, __GMBD_REAL(1), x, M, __GMBD_REAL(1));
}
/**
 * v = v+M*x
 * Will not resize v, as this would not make sense
 **/
inline static
void addmult(tVarVector<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  addmult(v, __GMBD_REAL(1), x, M, __GMBD_REAL(1));
}

/**
 * v = a+(x^T*M)^T
 **/
inline static
void addmult(tVectorView<__GMBD_REAL> &v,
             const tConstVectorView<__GMBD_REAL> &a, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  v = a;
  addmult(v, x, M);
}

/**
 * v = a+(x^T*M)^T
 **/
inline static
void addmult(tVarVector<__GMBD_REAL> &v,
             const tConstVectorView<__GMBD_REAL> &a, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  v = a;
  addmult(v, x, M);
}

/**
 * v = alpha*x^T*M
 **/
inline static
void mult(tVectorView<__GMBD_REAL> &v,
          const __GMBD_REAL& alpha, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  ::gravis::matrix::clear(v);
  addmult(v, alpha, x, M, 1);
}
/**
 * v = alpha*x^T*M
 * Will not resize v, as this would not make sense
 **/
inline static
void mult(tVarVector<__GMBD_REAL> &v,
          const __GMBD_REAL& alpha, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  ::gravis::matrix::clear(v);
  addmult(v, alpha, x, M, 1);
}

/**
 * v = x^T*M
 **/
inline static
void mult(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  mult(v, 1, x, M);
}

/**
 * v = x^T*M
 **/
inline static
void mult(tVarVector<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  mult(v, 1, x, M);
}

static inline
__GMBD_REAL abs(const __GMBD_REAL& a)
{
  return a< __GMBD_REAL(0) ? -a : a;
}

//// Norms
/** l1 norm **/
inline static __GMBD_REAL normL1(const tConstVectorView<__GMBD_REAL> &v)
{
  if (v.size()==0) return 0;
  __GMBD_REAL result = abs(v[0]);
  for (size_t i=1; i<v.size(); ++i) result += abs(v[i]);
  return result;
}
/** l1 norm **/
inline static __GMBD_REAL normL1(const tConstMatrixView<__GMBD_REAL> &v)
{
  if (v.size()==0) return 0;
  __GMBD_REAL result = abs(v[0]);
  for (size_t i=1; i<v.size(); ++i) result += abs(v[i]);
  return result;
}

/** l2 norm **/
inline static __GMBD_REAL normL2(const tConstVectorView<__GMBD_REAL> &v)
{
  return v.size()==0 ? 0 : __GMBD_xNRM2(v.size(), v.data, 1);
}
/** Frobenius norm **/
inline static __GMBD_REAL normL2(const tConstMatrixView<__GMBD_REAL> &v)
{
  return v.size()==0 ? 0 : __GMBD_xNRM2(v.size(), v.data, 1);
}

/** Squared l2 norm **/
inline static __GMBD_REAL normL2sqr(const tConstVectorView<__GMBD_REAL> &v)
{
  return ::gravis::matrix::priv::sqr(normL2(v));
}
/** Squared Frobenius norm **/
inline static __GMBD_REAL normL2sqr(const tConstMatrixView<__GMBD_REAL> &v)
{
  return ::gravis::matrix::priv::sqr(normL2(v));
}

/** linf norm **/
inline static __GMBD_REAL normLinf(const tConstVectorView<__GMBD_REAL> &v)
{
  if (v.size()==0) return 0;
  __GMBD_REAL result = abs(v[0]);
  for (size_t i=1; i<v.size(); ++i) result = std::max(result, abs(v[i]));
  return result;
}
/** linf norm **/
inline static __GMBD_REAL normLinf(const tConstMatrixView<__GMBD_REAL> &v)
{
  if (v.size()==0) return 0;
  __GMBD_REAL result = abs(v[0]);
  for (size_t i=1; i<v.size(); ++i) result = std::max(result, abs(v[i]));
  return result;
}

//// Matrix Matrix Operations
/** v += u **/
inline static void add(tVectorView<__GMBD_REAL> &v, const __GMBD_REAL& s, const tConstVectorView<__GMBD_REAL> &u)
{
  __GMBD_xAXPY(v.size(), s, u.data, 1, v.data, 1);
}
/** v += u **/
inline static void add(tVarVector<__GMBD_REAL>  &v, const __GMBD_REAL& s, const tConstVectorView<__GMBD_REAL> &u)
{
  __GMBD_xAXPY(v.size(), s, u.data, 1, v.data, 1);
}

/** V += U **/
inline static void add(tMatrixView<__GMBD_REAL> &V, const __GMBD_REAL& s, const tConstMatrixView<__GMBD_REAL> &U)
{
  __GMBD_xAXPY(V.size(), s, U.data, 1, V.data, 1);
}
/** V += U **/
inline static void add(tVarMatrix<__GMBD_REAL>  &V, const __GMBD_REAL& s, const tConstMatrixView<__GMBD_REAL> &U)
{
  __GMBD_xAXPY(V.size(), s, U.data, 1, V.data, 1);
}

/** v += u **/
inline static void add(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &u)
{
  add(v, 1, u);
}
/** v += u **/
inline static void add(tVarVector<__GMBD_REAL>  &v, const tConstVectorView<__GMBD_REAL> &u)
{
  add(v, 1, u);
}

/** V += U **/
inline static void add(tMatrixView<__GMBD_REAL> &V, const tConstMatrixView<__GMBD_REAL> &U)
{
  add(V, 1, U);
}
/** V += U **/
inline static void add(tVarMatrix<__GMBD_REAL>  &V, const tConstMatrixView<__GMBD_REAL> &U)
{
  add(V, 1, U);
}

/** r = v + u **/
inline static void add(tVectorView<__GMBD_REAL> &r, const tConstVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &u)
{
  r = v;
  add(r,u);
}
/** r = v + u **/
inline static void add(tVarVector<__GMBD_REAL>  &r, const tConstVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &u)
{
  r = v;
  add(r,u);
}

/** R = V + U **/
inline static void add(tMatrixView<__GMBD_REAL> &R, const tConstMatrixView<__GMBD_REAL> &V, const tConstMatrixView<__GMBD_REAL> &U)
{
  R = V;
  add(R,U);
}
/** R = V + U **/
inline static void add(tVarMatrix<__GMBD_REAL>  &R, const tConstMatrixView<__GMBD_REAL> &V, const tConstMatrixView<__GMBD_REAL> &U)
{
  R = V;
  add(R,U);
}


/** v -= u **/
inline static void sub(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &u)
{
  add(v, -1, u);
}
/** v -= u **/
inline static void sub(tVarVector<__GMBD_REAL>  &v, const tConstVectorView<__GMBD_REAL> &u)
{
  add(v, -1, u);
}

/** V -= U **/
inline static void sub(tMatrixView<__GMBD_REAL> &V, const tConstMatrixView<__GMBD_REAL> &U)
{
  add(V, -1, U);
}
/** V -= U **/
inline static void sub(tVarMatrix<__GMBD_REAL>  &V, const tConstMatrixView<__GMBD_REAL> &U)
{
  add(V, -1, U);
}

/** v -= su **/
inline static void sub(tVectorView<__GMBD_REAL> &v, const __GMBD_REAL& s, const tConstVectorView<__GMBD_REAL> &u)
{
  add(v, -s, u);
}
/** v -= su **/
inline static void sub(tVarVector<__GMBD_REAL>  &v, const __GMBD_REAL& s, const tConstVectorView<__GMBD_REAL> &u)
{
  add(v, -s, u);
}

/** V -= sU **/
inline static void sub(tMatrixView<__GMBD_REAL> &V, const __GMBD_REAL& s, const tConstMatrixView<__GMBD_REAL> &U)
{
  add(V, -s, U);
}
/** V -= sU **/
inline static void sub(tVarMatrix<__GMBD_REAL>  &V, const __GMBD_REAL& s, const tConstMatrixView<__GMBD_REAL> &U)
{
  add(V, -s, U);
}

//// Matrix Scalar Operations
/** Arithmethic operations with scalars *= **/
inline static void mult(tVectorView<__GMBD_REAL> &o, const tConstVectorView<__GMBD_REAL> &v, const __GMBD_REAL& scalar)
{
  o=v;
  __GMBD_xSCAL(v.size(), scalar, o.data, 1);
}
/** Arithmethic operations with scalars *= **/
inline static void mult(tVarVector<__GMBD_REAL>  &o, const tConstVectorView<__GMBD_REAL> &v, const __GMBD_REAL& scalar)
{
  o=v;
  __GMBD_xSCAL(v.size(), scalar, o.data, 1);
}

/** Arithmethic operations with scalars *= **/
inline static void mult(tMatrixView<__GMBD_REAL> &o, const tConstMatrixView<__GMBD_REAL> &v, const __GMBD_REAL& scalar)
{
  o=v;
  __GMBD_xSCAL(v.size(), scalar, o.data, 1);
}
/** Arithmethic operations with scalars *= **/
inline static void mult(tVarMatrix<__GMBD_REAL>  &o, const tConstMatrixView<__GMBD_REAL> &v, const __GMBD_REAL& scalar)
{
  o=v;
  __GMBD_xSCAL(v.size(), scalar, o.data, 1);
}

/** Arithmethic operations with scalars = * **/
inline static void mult(tVectorView<__GMBD_REAL> &v, const __GMBD_REAL& scalar)
{
  __GMBD_xSCAL(v.size(), scalar, v.data, 1);
}
/** Arithmethic operations with scalars = * **/
inline static void mult(tVarVector<__GMBD_REAL>  &v, const __GMBD_REAL& scalar)
{
  __GMBD_xSCAL(v.size(), scalar, v.data, 1);
}

/** Arithmethic operations with scalars = * **/
inline static void mult(tMatrixView<__GMBD_REAL> &v, const __GMBD_REAL& scalar)
{
  __GMBD_xSCAL(v.size(), scalar, v.data, 1);
}
/** Arithmethic operations with scalars = * **/
inline static void mult(tVarMatrix<__GMBD_REAL>  &v, const __GMBD_REAL& scalar)
{
  __GMBD_xSCAL(v.size(), scalar, v.data, 1);
}

/** Dotproduct * **/
inline static __GMBD_REAL dot(const tVarVector<__GMBD_REAL>  &u, const tVarVector<__GMBD_REAL>  &v)
{
  return __GMBD_xDOT(u.size(),u.data,1,v.data,1);
}
inline static __GMBD_REAL dot(const tVectorView<__GMBD_REAL>  &u, const tVectorView<__GMBD_REAL>  &v)
{
  return __GMBD_xDOT(u.size(),u.data,1,v.data,1);
}

inline static void pinv(tVarMatrix<__GMBD_REAL>  &A)
{
  tVarMatrix<__GMBD_REAL> U(A.h,A.h);
  tVarVector<__GMBD_REAL> s(std::min(A.w,A.h));
  tVarMatrix<__GMBD_REAL> VT(A.w,A.w);
  svd_dc(U, s, VT, A);

  for (unsigned int i = 0; i < s.h; ++i)
  {
    if(s[i] != 0) s[i] = 1 / s[i];
  }
  //GEMM(VT,A)
  // UT = diag(s) * UT
  for (unsigned int i = 0; i < U.w; ++i)
  {
    for (unsigned int j = 0; j < U.h; ++j)
    {
      U(j,i) = (i < s.h) ? s(i) * U(j,i) : 0;
    }
  }

  // A = V * UT
  //tVarMatrix<__GMBD_REAL> X(VT.w, U.h);
  A.resize(A.w,A.h);
  for (unsigned int i = 0; i < A.h; ++i)
  {
    for (unsigned int j = 0; j < A.w; ++j)
    {
      A(i,j) = 0;
      for (unsigned int k = 0; k < std::min(VT.h,U.h); ++k)
      {
        A(i,j) += VT(k,i) * U(j,k);
      }
    }
  }
  //A = X;
}

#undef __GMBD_REAL
#undef __GMBD_xGEMV
#undef __GMBD_xNRM2
#undef __GMBD_xSCAL
#undef __GMBD_xAXPY
#undef __GMBD_xDOT
#undef __GMBD_xGESVD
#undef __GMBD_xGESDD
