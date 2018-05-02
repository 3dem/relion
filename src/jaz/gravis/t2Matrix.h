#ifndef __LIBGRAVIS_T2MATRIX_H__
#define __LIBGRAVIS_T2MATRIX_H__
/******************************************************************************
**        Title: t2Matrix.h
**  Description: Represents a 2x2 matrix with column-major memory layout.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Brian Amberg 2006
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <iostream>
#include "t2Vector.h"

namespace gravis
{

  /**
   * A 2x2 matrix with column-major memory layout
   **/
  template <class T>
  struct t2Matrix
  {

    T m[4];
    t2Matrix ()
    {
      loadIdentity();
    };
    t2Matrix (const t2Matrix& mat);
    explicit t2Matrix (const T& val)
    {
      m[0] = val;
      m[1] = val;
      m[2] = val;
      m[3] = val;
    };
    explicit t2Matrix (const T* v_ptr)
    {
      m[ 0] = v_ptr[ 0];
      m[ 1] = v_ptr[ 1];
      m[ 2] = v_ptr[ 2];
      m[ 3] = v_ptr[ 3];
    };

    t2Matrix (T m0, T m2,
              T m1, T m3);

    void set (T m0, T m2,
              T m1, T m3);

    const T& operator[] (int idx) const;
    T& operator[] (int idx);

    const T& operator() (int row, int col) const
    {
      return m[(col << 1) + row];
    }
    T& operator() (int row, int col)
    {
      return m[(col << 1) + row];
    }

    t2Matrix operator*(T f) const;

    t2Vector<T> operator*  (const t2Vector<T>&) const;

    t2Matrix    operator*  (const t2Matrix&) const;
    t2Matrix&   operator*= (const t2Matrix&);

    /**
     * Element Wise Addition (Inplace)
     **/
    inline t2Matrix&   operator+= (const t2Matrix& rhs)
    {
      for (size_t i=0; i<4; ++i) m[i] += rhs.m[i];
      return *this;
    }
    /**
     * Element Wise Subtraction (Inplace)
     **/
    inline t2Matrix&   operator-= (const t2Matrix& rhs)
    {
      for (size_t i=0; i<4; ++i) m[i] -= rhs.m[i];
      return *this;
    }
    /**
     * Element Wise Addition
     **/
    inline t2Matrix    operator+  (const t2Matrix& rhs) const
    {
      t2Matrix result(*this);
      return(result += rhs);
    }
    /**
     * Element Wise Subtraction
     **/
    inline t2Matrix    operator-  (const t2Matrix& rhs) const
    {
      t2Matrix result(*this);
      return(result -= rhs);
    }

    /**
     * Matrix Norm (2 Norm)
     **/
    inline T norm2() const
    {
      return m[0]*m[0] + m[1]*m[1] + m[2]*m[2] + m[3]*m[3];
    }

    /**
     * Matrix Trace (sum(diag(M)))
     **/
    T trace() const;

    void transpose ();
    void invert();
    void loadIdentity ();

    static t2Matrix identity();

    static t2Matrix scale (const t2Vector<T>& s)
    {
      return t2Matrix<T>(
               s[0],   T(0),
               T(0)  , s[1]);
    };
    static t2Matrix scale (const T& s)
    {
      return scale(t2Vector<T>(s,s));
    }

    /**
     * Create a 2x2 rotation matrix for a clockwise rotation around a rad.
     **/
    static t2Matrix rotation(const T& a)
    {
      return t2Matrix(
               cos(a), -sin(a),
               sin(a), cos(a));
    }



  };

  template <class T> inline
  t2Matrix<T>::t2Matrix (const t2Matrix& mat)
  {
    m[0] = mat.m[0];
    m[1] = mat.m[1];
    m[2] = mat.m[2];
    m[3] = mat.m[3];
  }


  template <class T> inline
  t2Matrix<T>::t2Matrix (T m0, T m2,
                         T m1, T m3)
  {
    m[ 0] =  m0;
    m[ 1] =  m1;
    m[ 2] =  m2;
    m[ 3] =  m3;
  }

  template <class T> inline
  void t2Matrix<T>::set (T m0, T m2,
                         T m1, T m3)
  {
    m[ 0] =  m0;
    m[ 1] =  m1;
    m[ 2] =  m2;
    m[ 3] =  m3;
  }

  template <class T> inline
  const T& t2Matrix<T>::operator[] (int idx) const
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert((idx >= 0) && (idx < 4));
#endif
    return m[idx];
  }


  template <class T> inline
  T& t2Matrix<T>::operator[] (int idx)
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert((idx >= 0) && (idx < 4));
#endif
    return m[idx];
  }


  template <class T> inline
  t2Matrix<T>& t2Matrix<T>::operator*= (const t2Matrix<T>& op)
  {
    *this = this->operator*(op);
    return *this;
  }


  template <class T> inline
  t2Matrix<T> t2Matrix<T>::operator* (const t2Matrix<T>& op) const
  {
    return t2Matrix(
             m[0]*op.m[ 0] + m[2]*op.m[ 1], m[0]*op.m[ 2] + m[2]*op.m[3],
             m[1]*op.m[ 0] + m[3]*op.m[ 1], m[1]*op.m[ 2] + m[3]*op.m[3]);
  }


  template <class T> inline
  t2Vector<T> t2Matrix<T>::operator* (const t2Vector<T>& op) const
  {
    return t2Vector<T>(
             m[ 0]*op.x + m[ 2]*op.y, m[ 1]*op.x + m[ 3]*op.y);
  }

  template <class T>
  inline
  t2Matrix<T> t2Matrix<T>::operator* (T f) const
  {
    return t2Matrix<T>(f * m[0], f * m[2], f * m[1], f * m[3]);
  }

  template <class T>
  inline
  t2Matrix<T> operator* (T f, const t2Matrix<T>& v)
  {
    return t2Matrix<T>(f * v[0], f * v[2], f * v[1], f * v[3]);
  }


  template <class T> inline
  void t2Matrix<T>::loadIdentity ()
  {
    m[ 0] = T(1);
    m[ 2] = T(0);
    m[ 1] = T(0);
    m[ 3] = T(1);
  }

  template <class T> inline
  void t2Matrix<T>::transpose ()
  {
    std::swap(m[1], m[2]);
  }

  template <class T> inline
  void t2Matrix<T>::invert()
  {
    t2Matrix<T> A = *this;
    T di = 1.0/(A[0]*A[3]-A[1]*A[2]);
    m[0]= A[3]*di;
    m[1]=-A[1]*di;
    m[2]=-A[2]*di;
    m[3]= A[0]*di;
  }

  template <class T> inline
  t2Matrix<T> t2Matrix<T>::identity ()
  {
    return t2Matrix<T>(
             T(1), T(0),
             T(0), T(1));
  }

  template <class T> inline
  T t2Matrix<T>::trace() const
  {
    return ( m[0] + m[3]);
  }


  template <class T> inline
  std::ostream& operator<< (std::ostream& os, const t2Matrix<T>& arg)
  {
    os << "[ " << arg[ 0] << "  " << arg[ 2]  << " ]\n";
    os << "[ " << arg[ 1] << "  " << arg[ 3]  << " ]\n";
    return os;
  }

  template <class T> inline
  std::istream& operator>> ( std::istream& is, t2Matrix<T>& arg)
  {
    std::string dummy;
    is >> dummy >> arg[ 0] >> arg[ 2] >> dummy;
    is >> dummy >> arg[ 1] >> arg[ 3] >> dummy;
    return is;
  }

  typedef t2Matrix<float > f2Matrix;
  typedef t2Matrix<double> d2Matrix;

} /* Close Namespace "gravis" */

#endif
