#ifndef __LIBGRAVIS_T4MATRIX_H__
#define __LIBGRAVIS_T4MATRIX_H__
/******************************************************************************
 **        Title: t4Matrix.h
 **  Description: Represents a 4x4 matrix with column-major memory layout.
 **
 **       Author: Jean-Sebastien Pierrard, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#include <iostream>
#include "t4Vector.h"
#include "private/tDeterminants.h"
#include "t3Matrix.h"

namespace gravis
{

  /*! \brief A 4x4 matrix class.
   *
   * There is no operator*=, because some people expect it to be a left-multiplication
   * and others a right-multiplication. To avoid confusion we only provide the explicit
   * methods lmul() and rmul().
   */
  template <class T>
  class t4Matrix
  {

    public:
      T m[16];

      typedef T scalar_type;

      t4Matrix();
      explicit t4Matrix(T val);
      t4Matrix(const T* v_ptr);
      t4Matrix(const t4Matrix& mat);
      t4Matrix(const t3Matrix<T>& mat);

      t4Matrix(T m0, T m4, T  m8, T m12,
               T m1, T m5, T  m9, T m13,
               T m2, T m6, T m10, T m14,
               T m3, T m7, T m11, T m15);

      template <class S>
      explicit t4Matrix (const t4Matrix<S>& mat)
      {
        for(int i=0; i<16; ++i)
          m[i] = static_cast<T>(mat.m[i]);
      }

      void set( T m0, T m4, T  m8, T m12,
                T m1, T m5, T  m9, T m13,
                T m2, T m6, T m10, T m14,
                T m3, T m7, T m11, T m15);

      //! Check if the entries of the other vector differ by less than epsilon.
      //  It is better to use this than to use operator== for comparision, if it is
      //  not the same vertex.
      bool isClose( const t4Matrix& o, const T epsilon) const
      {
        for (int i=0; i<16; i++)
          if (std::fabs(m[i]-o.m[i]) >= epsilon)
            return false;
        return true;
      }


      bool operator==(const t4Matrix& o) const
      {
        for (int i=0; i<16; i++)
          if (m[i] != o.m[i]) return false;
        return true;
      }

      bool operator!=(const t3Matrix<T> &o) const
      {
        return !(*this == o);
      }

      /*! \brief Return indexed entry (column major). */
      const T& operator[] (int idx) const
      {
        return m[idx];
      }
      /*! \brief Return reference to indexed entry (column major). */
      T& operator[] (int idx)
      {
        return m[idx];
      }

      /*! \brief Return entry in row i and column j. */
      const T& operator() (int row, int col) const
      {
        return m[col * 4 + row];
      }
      /*! \brief Return reference to entry in row i and column j. */
      T& operator() (int row, int col)
      {
        return m[col * 4 + row];
      }

      t4Matrix operator*(T f) const;
      t4Matrix& operator*=(T f);
      t4Matrix operator/(T f) const;
      t4Matrix& operator/=(T f);
      t4Vector<T> operator*(const t4Vector<T>&) const;
      t4Matrix operator*(const t4Matrix&) const;
      t4Matrix& operator+=(const t4Matrix&);
      t4Matrix& operator-=(const t4Matrix&);
      t4Matrix operator+(const t4Matrix&) const;
      t4Matrix operator-(const t4Matrix&) const;
      t4Matrix operator-() const;

      t4Matrix& lmul(const t4Matrix& m);
      t4Matrix& rmul(const t4Matrix& m);

      T trace() const;
	  T det() const;
	  T FrobeniusNorm() const;

	  t4Matrix adjugate() const;
      t4Matrix& transpose();
      t4Matrix& invert();
      t4Matrix& loadIdentity();
      t4Matrix& copy(const t3Matrix<T>& mat);

      static t4Matrix translation(const t3Vector<T>&);
      static t4Matrix scale(const t3Vector<T>&);
      static t4Matrix scale (const T& s)
      {
        return scale(t3Vector<T>(s,s,s));
      }
      static t4Matrix rotation(const t3Vector<T>& u, const t3Vector<T>& v);
      static t4Matrix rotation(const t3Vector<T>& axis, float angle);
      static t4Matrix rotationX(T angle);
      static t4Matrix rotationY(T angle);
      static t4Matrix rotationZ(T angle);
  };

  /*! \brief Constructs an identity matrix. */
  template <class T> inline t4Matrix<T>::t4Matrix()
  {
    loadIdentity();
  }

  /*! \brief Constructs a matrix with all entries set to val. */
  template <class T> inline
  t4Matrix<T>::t4Matrix(T val)
  {
    for (int i = 0; i < 16; i++) m[i] = val;
  }

  /*! \brief Constructs a matrix with entries taken from an array.
   *
   * \param v_ptr array must be of appropriate length and in column-major layout */
  template <class T> inline
  t4Matrix<T>::t4Matrix(const T* v_ptr)
  {
    for (int i = 0; i < 16; i++) m[i] = v_ptr[i];
  }

  /*! \brief Copy constructor. */
  template <class T> inline
  t4Matrix<T>::t4Matrix(const t4Matrix& mat)
  {
    for (int i = 0; i < 16; i++) m[i] = mat.m[i];
  }

  /*! \brief Copy constructor. */
  template <class T> inline
  t4Matrix<T>::t4Matrix(const t3Matrix<T>& mat)
  {

    m[ 0] = mat[0];
    m[ 1] = mat[1];
    m[ 2] = mat[2];
    m[ 3] = 0.f;
    m[ 4] = mat[3];
    m[ 5] = mat[4];
    m[ 6] = mat[5];
    m[ 7] = 0.f;
    m[ 8] = mat[6];
    m[ 9] = mat[7];
    m[10] = mat[8];
    m[11] = 0.f;
    m[12] = 0.f;
    m[13] = 0.f;
    m[14] = 0.f;
    m[15] = 1.f;
  }

  /*! \brief Constructs a matrix from the given entries (row major). */
  template <class T> inline
  t4Matrix<T>::t4Matrix (T m0, T m4, T  m8, T m12,
                         T m1, T m5, T  m9, T m13,
                         T m2, T m6, T m10, T m14,
                         T m3, T m7, T m11, T m15)
  {
    m[ 0] =  m0;
    m[ 1] =  m1;
    m[ 2] =  m2;
    m[ 3] =  m3;
    m[ 4] =  m4;
    m[ 5] =  m5;
    m[ 6] =  m6;
    m[ 7] =  m7;
    m[ 8] =  m8;
    m[ 9] =  m9;
    m[10] = m10;
    m[11] = m11;
    m[12] = m12;
    m[13] = m13;
    m[14] = m14;
    m[15] = m15;
  }

  /*! \brief Overwrites this matrix with the given entries (row major). */
  template <class T> inline
  void t4Matrix<T>::set (T m0, T m4, T  m8, T m12,
                         T m1, T m5, T  m9, T m13,
                         T m2, T m6, T m10, T m14,
                         T m3, T m7, T m11, T m15 )
  {
    m[ 0] =  m0;
    m[ 1] =  m1;
    m[ 2] =  m2;
    m[ 3] =  m3;
    m[ 4] =  m4;
    m[ 5] =  m5;
    m[ 6] =  m6;
    m[ 7] =  m7;
    m[ 8] =  m8;
    m[ 9] =  m9;
    m[10] = m10;
    m[11] = m11;
    m[12] = m12;
    m[13] = m13;
    m[14] = m14;
    m[15] = m15;
  }


  /*! \brief Scalar times matrix. */
  template <class T> inline
  t4Matrix<T> operator*(T f, const t4Matrix<T>& mat)
  {
    t4Matrix<T> out(mat);
    out *= f;
    return out;
  }

  /*! \brief Matrix times scalar. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::operator*(T f) const
  {
    t4Matrix<T> out(*this);
    out *= f;
    return out;
  }

  /*! \brief Multiply this matrix with a scalar. */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::operator*=(T f)
  {
    for (int i = 0; i < 16; i++) m[i] *= f;
    return *this;
  }

  /*! \brief Matrix divided by scalar. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::operator/(const T f) const
  {
    t4Matrix<T> out(*this);
    out /= f;
    return out;
  }

  /*! \brief Divide this matrix by a scalar. */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::operator/=(const T f)
  {
    for (int i = 0; i < 16; i++) m[i] /= f;
    return *this;
  }


  /*! \brief Matrix times vector. */
  template <class T> inline
  t4Vector<T> t4Matrix<T>::operator*(const t4Vector<T>& op) const
  {
    return t4Vector<T>(
             m[ 0]*op.x + m[ 4]*op.y + m[ 8]*op.z + m[12]*op.w,
             m[ 1]*op.x + m[ 5]*op.y + m[ 9]*op.z + m[13]*op.w,
             m[ 2]*op.x + m[ 6]*op.y + m[10]*op.z + m[14]*op.w,
             m[ 3]*op.x + m[ 7]*op.y + m[11]*op.z + m[15]*op.w
           );
  }

  /*! \brief Matrix times matrix. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::operator* (const t4Matrix<T>& op) const
  {

    return t4Matrix(
             m[0]*op.m[ 0] + m[4]*op.m[ 1] + m[8]*op.m[ 2] + m[12]*op.m[ 3],      // ROW 1
             m[0]*op.m[ 4] + m[4]*op.m[ 5] + m[8]*op.m[ 6] + m[12]*op.m[ 7],
             m[0]*op.m[ 8] + m[4]*op.m[ 9] + m[8]*op.m[10] + m[12]*op.m[11],
             m[0]*op.m[12] + m[4]*op.m[13] + m[8]*op.m[14] + m[12]*op.m[15],

             m[1]*op.m[ 0] + m[5]*op.m[ 1] + m[9]*op.m[ 2] + m[13]*op.m[ 3],      // ROW 2
             m[1]*op.m[ 4] + m[5]*op.m[ 5] + m[9]*op.m[ 6] + m[13]*op.m[ 7],
             m[1]*op.m[ 8] + m[5]*op.m[ 9] + m[9]*op.m[10] + m[13]*op.m[11],
             m[1]*op.m[12] + m[5]*op.m[13] + m[9]*op.m[14] + m[13]*op.m[15],

             m[2]*op.m[ 0] + m[6]*op.m[ 1] + m[10]*op.m[ 2] + m[14]*op.m[ 3],     // ROW 3
             m[2]*op.m[ 4] + m[6]*op.m[ 5] + m[10]*op.m[ 6] + m[14]*op.m[ 7],
             m[2]*op.m[ 8] + m[6]*op.m[ 9] + m[10]*op.m[10] + m[14]*op.m[11],
             m[2]*op.m[12] + m[6]*op.m[13] + m[10]*op.m[14] + m[14]*op.m[15],

             m[3]*op.m[ 0] + m[7]*op.m[ 1] + m[11]*op.m[ 2] + m[15]*op.m[ 3],     // ROW 4
             m[3]*op.m[ 4] + m[7]*op.m[ 5] + m[11]*op.m[ 6] + m[15]*op.m[ 7],
             m[3]*op.m[ 8] + m[7]*op.m[ 9] + m[11]*op.m[10] + m[15]*op.m[11],
             m[3]*op.m[12] + m[7]*op.m[13] + m[11]*op.m[14] + m[15]*op.m[15]
           );
  }

  /*! \brief Adds other matrix to this matrix. */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::operator+=(const t4Matrix<T>& op)
  {
    for (int i = 0; i < 16; i++) m[i] += op.m[i];
    return *this;
  }

  /*! \brief Subtracts other matrix from this matrix. */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::operator-=(const t4Matrix<T>& op)
  {
    *this += -op;
    return *this;
  }

  /*! \brief Matrix plus matrix. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::operator+(const t4Matrix<T>& op) const
  {
    t4Matrix<T> out(*this);
    return out += op;
  }

  /*! \brief Matrix minus matrix. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::operator-(const t4Matrix<T>& op) const
  {
    t4Matrix<T> out(*this);
    return out += -op;
  }

  /*! \brief Return additive inverse of this matrix. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::operator-() const
  {
    t4Matrix<T> out(*this);
    for (int i = 0; i < 16; i++) out[i] = -out[i];
    return out;
  }

  /*! \brief Right-multiply m to this matrix (*this = *this * m). */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::rmul(const t4Matrix<T>& m)
  {
    *this = *this * m;
    return *this;
  }

  /*! \brief Left-multiply m to this matrix (*this = m * *this). */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::lmul(const t4Matrix<T>& m)
  {
    *this = m * *this;
    return *this;
  }

  /*! \brief Return the trace of this matrix (\f$a_{11} + a_{22} + a_{33} + a_{44}\f$). */
  template <class T> inline
  T t4Matrix<T>::trace() const
  {
    return ( m[0] + m[5] + m[10] + m[15] );
  }
  
  /*! \brief Return the determinant of this matrix. */
  template <class T> inline
  T t4Matrix<T>::det() const
  {
    return det4x4(m[ 0], m[ 4], m[ 8], m[12],
                  m[ 1], m[ 5], m[ 9], m[13],
                  m[ 2], m[ 6], m[10], m[14],
                  m[ 3], m[ 7], m[11], m[15]);
  }
  
  template <class T> inline
  T t4Matrix<T>::FrobeniusNorm() const
  {
    double sum(0.0);
	
    for (int i = 0; i < 16; i++)
    {
      sum += m[i] * m[i];
    }
	
    return sum;
  }
  
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::adjugate() const
  {
    const t4Matrix<T>& A = *this;
	t4Matrix<T> B;

    B(0,0) =  det3x3(A(1,1), A(2,1), A(3,1), A(1,2), A(2,2), A(3,2), A(1,3), A(2,3), A(3,3));
    B(1,0) = -det3x3(A(1,0), A(2,0), A(3,0), A(1,2), A(2,2), A(3,2), A(1,3), A(2,3), A(3,3));
    B(2,0) =  det3x3(A(1,0), A(2,0), A(3,0), A(1,1), A(2,1), A(3,1), A(1,3), A(2,3), A(3,3));
    B(3,0) = -det3x3(A(1,0), A(2,0), A(3,0), A(1,1), A(2,1), A(3,1), A(1,2), A(2,2), A(3,2));

    B(0,1) = -det3x3(A(0,1), A(2,1), A(3,1), A(0,2), A(2,2), A(3,2), A(0,3), A(2,3), A(3,3));
    B(1,1) =  det3x3(A(0,0), A(2,0), A(3,0), A(0,2), A(2,2), A(3,2), A(0,3), A(2,3), A(3,3));
    B(2,1) = -det3x3(A(0,0), A(2,0), A(3,0), A(0,1), A(2,1), A(3,1), A(0,3), A(2,3), A(3,3));
    B(3,1) =  det3x3(A(0,0), A(2,0), A(3,0), A(0,1), A(2,1), A(3,1), A(0,2), A(2,2), A(3,2));

    B(0,2) =  det3x3(A(0,1), A(1,1), A(3,1), A(0,2), A(1,2), A(3,2), A(0,3), A(1,3), A(3,3));
    B(1,2) = -det3x3(A(0,0), A(1,0), A(3,0), A(0,2), A(1,2), A(3,2), A(0,3), A(1,3), A(3,3));
    B(2,2) =  det3x3(A(0,0), A(1,0), A(3,0), A(0,1), A(1,1), A(3,1), A(0,3), A(1,3), A(3,3));
    B(3,2) = -det3x3(A(0,0), A(1,0), A(3,0), A(0,1), A(1,1), A(3,1), A(0,2), A(1,2), A(3,2));

    B(0,3) = -det3x3(A(0,1), A(1,1), A(2,1), A(0,2), A(1,2), A(2,2), A(0,3), A(1,3), A(2,3));
    B(1,3) =  det3x3(A(0,0), A(1,0), A(2,0), A(0,2), A(1,2), A(2,2), A(0,3), A(1,3), A(2,3));
    B(2,3) = -det3x3(A(0,0), A(1,0), A(2,0), A(0,1), A(1,1), A(2,1), A(0,3), A(1,3), A(2,3));
    B(3,3) =  det3x3(A(0,0), A(1,0), A(2,0), A(0,1), A(1,1), A(2,1), A(0,2), A(1,2), A(2,2));

    return B;
  }

  /*! \brief Transpose this matrix.
   * Attention: Although innocent looking this is an inplace operation
   **/
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::transpose ()
  {
    std::swap(m[1], m[4]);
    std::swap(m[2], m[8]);
    std::swap(m[3], m[12]);
    std::swap(m[6], m[9]);
    std::swap(m[7], m[13]);
    std::swap(m[11], m[14]);
    return *this;
  }

  /*! \brief Invert this matrix.
   * Attention: Although innocent looking, this is an inplace operation
   **/
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::invert()
  {
    T det, oodet;
    t4Matrix<T> A = *this;

    (*this)(0,0) =  det3x3(A(1,1), A(2,1), A(3,1), A(1,2), A(2,2), A(3,2), A(1,3), A(2,3), A(3,3));
    (*this)(1,0) = -det3x3(A(1,0), A(2,0), A(3,0), A(1,2), A(2,2), A(3,2), A(1,3), A(2,3), A(3,3));
    (*this)(2,0) =  det3x3(A(1,0), A(2,0), A(3,0), A(1,1), A(2,1), A(3,1), A(1,3), A(2,3), A(3,3));
    (*this)(3,0) = -det3x3(A(1,0), A(2,0), A(3,0), A(1,1), A(2,1), A(3,1), A(1,2), A(2,2), A(3,2));

    (*this)(0,1) = -det3x3(A(0,1), A(2,1), A(3,1), A(0,2), A(2,2), A(3,2), A(0,3), A(2,3), A(3,3));
    (*this)(1,1) =  det3x3(A(0,0), A(2,0), A(3,0), A(0,2), A(2,2), A(3,2), A(0,3), A(2,3), A(3,3));
    (*this)(2,1) = -det3x3(A(0,0), A(2,0), A(3,0), A(0,1), A(2,1), A(3,1), A(0,3), A(2,3), A(3,3));
    (*this)(3,1) =  det3x3(A(0,0), A(2,0), A(3,0), A(0,1), A(2,1), A(3,1), A(0,2), A(2,2), A(3,2));

    (*this)(0,2) =  det3x3(A(0,1), A(1,1), A(3,1), A(0,2), A(1,2), A(3,2), A(0,3), A(1,3), A(3,3));
    (*this)(1,2) = -det3x3(A(0,0), A(1,0), A(3,0), A(0,2), A(1,2), A(3,2), A(0,3), A(1,3), A(3,3));
    (*this)(2,2) =  det3x3(A(0,0), A(1,0), A(3,0), A(0,1), A(1,1), A(3,1), A(0,3), A(1,3), A(3,3));
    (*this)(3,2) = -det3x3(A(0,0), A(1,0), A(3,0), A(0,1), A(1,1), A(3,1), A(0,2), A(1,2), A(3,2));

    (*this)(0,3) = -det3x3(A(0,1), A(1,1), A(2,1), A(0,2), A(1,2), A(2,2), A(0,3), A(1,3), A(2,3));
    (*this)(1,3) =  det3x3(A(0,0), A(1,0), A(2,0), A(0,2), A(1,2), A(2,2), A(0,3), A(1,3), A(2,3));
    (*this)(2,3) = -det3x3(A(0,0), A(1,0), A(2,0), A(0,1), A(1,1), A(2,1), A(0,3), A(1,3), A(2,3));
    (*this)(3,3) =  det3x3(A(0,0), A(1,0), A(2,0), A(0,1), A(1,1), A(2,1), A(0,2), A(1,2), A(2,2));

    det = (A(0,0) * (*this)(0,0)) + (A(0,1) * (*this)(1,0)) + (A(0,2) * (*this)(2,0)) + (A(0,3) * (*this)(3,0));

    oodet = T(1) / det;

    *this *= oodet;

    return *this;
  }

  /*! \brief Overwrite this matrix with an identity matrix. */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::loadIdentity ()
  {
    m[ 0] = T(1);
    m[ 1] = T(0);
    m[ 2] = T(0);
    m[ 3] = T(0);
    m[ 4] = T(0);
    m[ 5] = T(1);
    m[ 6] = T(0);
    m[ 7] = T(0);
    m[ 8] = T(0);
    m[ 9] = T(0);
    m[10] = T(1);
    m[11] = T(0);
    m[12] = T(0);
    m[13] = T(0);
    m[14] = T(0);
    m[15] = T(1);
    return *this;
  }

  /*! \brief Copies the 3x3 matrix into the upper left corner of this
   * instance. */
  template <class T> inline
  t4Matrix<T>& t4Matrix<T>::copy(const t3Matrix<T>& mat)
  {
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        m[4 * j + i] = mat.m[3 * j + i];
    return *this;
  }

  /*! \brief Return a matrix representing a translation by t. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::translation(const t3Vector<T>& t)
  {
    return t4Matrix<T>(T(1), T(0), T(0), t.x,
                       T(0), T(1), T(0), t.y,
                       T(0), T(0), T(1), t.z,
                       T(0), T(0), T(0), T(1)
                      );
  }

  /*! \brief Return a matrix represnting a scaling by s. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::scale (const t3Vector<T>& s)
  {
    return t4Matrix<T>(
             s.x,  T(0), T(0), T(0),
             T(0), s.y,  T(0), T(0),
             T(0), T(0), s.z,  T(0),
             T(0), T(0), T(0), T(1)
           );
  }

  /*! Return a matrix that will rotate u into v. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::rotation(const t3Vector<T>& u, const t3Vector<T>& v)
  {
    t4Matrix<T> out;
    out.copy(t3Matrix<T>::rotation(u, v));
    return out;
  }

  /*! \brief Return a matrix that rotates by specified angle (in degrees) around specified axis. */
  template <class T> inline
  t4Matrix<T> t4Matrix<T>::rotation(const t3Vector<T>& axis, float angle)
  {
    t4Matrix<T> out;
    out.copy(t3Matrix<T>::rotation(axis, angle));
    return out;
  }

  template <class T> inline
  t4Matrix<T> t4Matrix<T>::rotationX (T a)
  {
    return t4Matrix<T>(
             T(1), T(0),      T(0),       T(0),
             T(0), T(cos(a)), T(-sin(a)), T(0),
             T(0), T(sin(a)), T(cos(a)),  T(0),
             T(0), T(0),      T(0),       T(1)
           );
  }

  template <class T> inline
  t4Matrix<T> t4Matrix<T>::rotationY (T a)
  { // ATTENTION!!! This is actually wrong!, -sin is in the first column
    // but this could disrupt everything!!!
    // Sandro Schoenborn, 2013-04-09, sandro.schoenborn@unibas.ch
    // Clemens Blumer, 2013-04-09, clemens.blumer@unibas.ch
    return t4Matrix<T>(
             T(cos(a)), T(0), T(-sin(a)), T(0),
             T(0),      T(1), T(0),       T(0),
             T(sin(a)), T(0), T(cos(a)),  T(0),
             T(0),      T(0), T(0),       T(1)
           );
  }

  template <class T> inline
  t4Matrix<T> t4Matrix<T>::rotationZ (T a)
  {
    return t4Matrix<T>(
             T(cos(a)), T(-sin(a)), T(0), T(0),
             T(sin(a)), T(cos(a)),  T(0), T(0),
             T(0),      T(0),       T(1), T(0),
             T(0),      T(0),       T(0), T(1)
           );
  }

  // TODO: Set Fixed Precision
  template <class T> inline
  std::ostream& operator<< (std::ostream& os, const t4Matrix<T>& arg)
  {
    os << "[ " << arg[ 0] << "  " << arg[ 4] << "  " << arg[ 8] << "  " << arg[12] << " ]\n";
    os << "| " << arg[ 1] << "  " << arg[ 5] << "  " << arg[ 9] << "  " << arg[13] << " |\n";
    os << "| " << arg[ 2] << "  " << arg[ 6] << "  " << arg[10] << "  " << arg[14] << " |\n";
    os << "[ " << arg[ 3] << "  " << arg[ 7] << "  " << arg[11] << "  " << arg[15] << " ]\n";
    return os;
  }

  template <class T> inline
  std::istream& operator>> ( std::istream& is, t4Matrix<T>& arg)
  {
    std::string dummy;
    is >> dummy >> arg[ 0] >> arg[ 4] >> arg[ 8] >> arg[12] >> dummy;
    is >> dummy >> arg[ 1] >> arg[ 5] >> arg[ 9] >> arg[13] >> dummy;
    is >> dummy >> arg[ 2] >> arg[ 6] >> arg[10] >> arg[14] >> dummy;
    is >> dummy >> arg[ 3] >> arg[ 7] >> arg[11] >> arg[15] >> dummy;
    return is;
  }

  typedef t4Matrix<float> f4Matrix;
  typedef t4Matrix<double> d4Matrix;

} // namespace gravis

#endif
