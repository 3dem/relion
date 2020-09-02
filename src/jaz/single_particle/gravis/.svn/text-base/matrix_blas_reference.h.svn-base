// v = v+M*x
inline static
void addmult(tVectorView<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  GRAVIS_CHECK( v.size()  == M.h, "v and M are incompatible");
  GRAVIS_CHECK( x.size()  == M.w, "M and x are incompatible");
  for (size_t j=0; j<M.w; ++j)
    for (size_t i=0; i<M.h; ++i)
      v[i] += x[j]*M(i,j);
}

// v = a+M*x
inline static
void addmult(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &a, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  GRAVIS_CHECK( v.size() == a.size(), "v and a are incompatible");
  // Addition
  v = a;
  addmult(v, M, x);
}

// v = M*x
inline static
void mult(tVectorView<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  ::gravis::matrix::clear(v);
  addmult(v, M, x);
}

// v = v+(x^T M)^T
inline static
void addmult(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  GRAVIS_CHECK( v.size()  == M.w, "v and M are incompatible");
  GRAVIS_CHECK( x.size()  == M.h, "M and x are incompatible");
  for (size_t i=0; i<M.h; ++i)
    for (size_t j=0; j<M.w; ++j)
      v[j] += x[i]*M(i,j);
}

// v = a+(x^T M)^T
inline static
void addmult(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &a, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  GRAVIS_CHECK( v.size() == a.size(), "v and a are incompatible");
  // Addition
  v = a;
  addmult(v, M, x);
}

// v = (x^T M)^T
inline static
void mult(tVectorView<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &x, const tConstMatrixView<__GMBD_REAL> &M)
{
  ::gravis::matrix::clear(v);
  addmult(v, M, x);
}

/**
 * Squared l2 norm
 **/
inline static
__GMBD_REAL normL2sqr(const tConstVectorView<__GMBD_REAL> &v)
{
  if (v.size() == 0) return 0;
  __GMBD_REAL result = v[0]*v[0];
  for (size_t i=1; i<v.size(); ++i) result += v[i]*v[i];
  return result;
}
/**
 * Squared frobenius norm
 **/
inline static
__GMBD_REAL normL2sqr(const tConstMatrixView<__GMBD_REAL> &v)
{
  if (v.size() == 0) return 0;
  __GMBD_REAL result = v[0]*v[0];
  for (size_t i=1; i<v.size(); ++i) result += v[i]*v[i];
  return result;
}

/** l2 norm **/
inline static __GMBD_REAL normL2(const tConstVectorView<__GMBD_REAL> &v)
{
  return __GMBD_REAL(sqrt(normL2sqr(v)));
}
/** forbenius norm **/
inline static __GMBD_REAL normL2(const tConstMatrixView<__GMBD_REAL> &v)
{
  return __GMBD_REAL(sqrt(normL2sqr(v)));
}

// v = v+M*x
inline static
void addmult(tVarVector<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  tVectorView<__GMBD_REAL> vv(v);
  addmult(vv, M, x);
}

// v = a+M*x
inline static
void addmult(tVarVector<__GMBD_REAL> &v, const tConstVectorView<__GMBD_REAL> &a, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  tVectorView<__GMBD_REAL> vv(v);
  addmult(vv, a, M, x);
}

// v = M*x
inline static
void mult(tVarVector<__GMBD_REAL> &v, const tConstMatrixView<__GMBD_REAL> &M, const tConstVectorView<__GMBD_REAL> &x)
{
  tVectorView<__GMBD_REAL> vv(v);
  addmult(vv, M, x);
}
