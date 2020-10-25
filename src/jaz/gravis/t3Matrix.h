#ifndef __LIBGRAVIS_T3MATRIX_H__
#define __LIBGRAVIS_T3MATRIX_H__
/******************************************************************************
 **        Title: t3Matrix.h
 **  Description: Represents a 3x3 matrix with column-major memory layout.
 **
 **       Author: Michael Keller, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#include <iostream>
#include <iomanip>
#include "t3Vector.h"
#include "tRGB.h"
#include "private/tDeterminants.h"

namespace gravis
{

  template <class T>
  class t4Matrix;

  /*! \brief A 3x3 matrix class.
   *
   * There is no operator*=, because some people expect it to be a left-multiplication
   * and others a right-multiplication. To avoid confusion we only provide the explicit
   * methods lmul() and rmul().
   */
  template <class T>
  class t3Matrix
  {

    public:
      T m[9];

      t3Matrix();
      explicit t3Matrix(T v);
      t3Matrix(const T* v_ptr);
      t3Matrix(const t3Matrix& mat);
      t3Matrix(T m0, T m3, T m6,
               T m1, T m4, T m7,
               T m2, T m5, T m8);

      template <class S>
      explicit t3Matrix (const t3Matrix<S>& mat)
      {
        for(int i=0; i<9; ++i)
          m[i] = static_cast<T>(mat.m[i]);
      }

      void set(T m0, T m3, T m6,
               T m1, T m4, T m7,
               T m2, T m5, T m8);

      /*! \brief Return indexed entry (column major). */
      const T& operator[](int idx) const
      {
        return m[idx];
      }
      /*! \brief Return reference to indexed entry (column major). */
      T& operator[](int idx)
      {
        return m[idx];
      }

      /*! \brief Return entry in row i and column j. */
      const T& operator()(int row, int col) const
      {
        return m[col * 3 + row];
      }
      /*! \brief Return reference to entry in row i and column j. */
      T& operator()(int row, int col)
      {
        return m[col * 3 + row];
      }


      //! Check if the entries of the other vector differ by less than epsilon.
      //  It is better to use this than to use operator== for comparision, if it is
      //  not the same vertex.
      bool isClose( const t3Matrix& o, const T epsilon) const
      {
        for (int i=0; i<9; i++)
          if (std::fabs(m[i]-o.m[i]) >= epsilon)
            return false;
        return true;
      }

      bool operator==(const t3Matrix& o) const
      {
        for (int i=0; i<9; i++)
          if (m[i] != o.m[i]) return false;
        return true;
      }

      bool operator!=(const t3Matrix& o) const
      {
        return !(*this == o);
      }

	  t3Matrix operator*(T f) const;
	  t3Matrix operator/(T f) const;
      t3Matrix& operator*=(const T& f);
      t3Matrix& operator/=(const T& f);
      t3Vector<T> operator*(const t3Vector<T>&) const;
      tRGB<T> operator*(const tRGB<T>&) const;
      t3Matrix operator*(const t3Matrix&) const;
      t3Matrix& operator+=(const t3Matrix&);
      t3Matrix& operator-=(const t3Matrix&);
      t3Matrix operator+(const t3Matrix&) const;
      t3Matrix operator-(const t3Matrix&) const;
      t3Matrix operator-() const;

      t3Matrix& lmul(const t3Matrix& m);
      t3Matrix& rmul(const t3Matrix& m);

      T trace() const;
      T det() const;

      t3Matrix adjugate() const;
      t3Matrix& transpose();
      t3Matrix& invert();
      t3Matrix& loadIdentity();

      t3Vector<T> getAxis() const;

      static t3Matrix extract(const t4Matrix<T>& mat);
      static t3Matrix scale(const t3Vector<T>&);
      static t3Matrix scale (const T& s)
      {
        return scale(t3Vector<T>(s,s,s));
      }
      static t3Matrix rotation(const t3Vector<T>& u, const t3Vector<T>& v);
      static t3Matrix rotation(const t3Vector<T>& axis, float angle);
      static t3Matrix rotationX(T angle);
      static t3Matrix rotationY(T angle);
      static t3Matrix rotationZ(T angle);
  };


  /*! \brief Constructs an identity matrix. */
  template <class T> inline t3Matrix<T>::t3Matrix()
  {
    loadIdentity();
  }

  /*! \brief Constructs a matrix with all entries set to val. */
  template <class T> inline
  t3Matrix<T>::t3Matrix(T val)
  {
    for (int i = 0; i < 9; i++) m[i] = val;
  }

  /*! \brief Constructs a matrix with entries taken from an array.
   *
   * \param v_ptr array must be of appropriate length and in column-major layout */
  template <class T> inline
  t3Matrix<T>::t3Matrix(const T* v_ptr)
  {
    for (int i = 0; i < 9; i++) m[i] = v_ptr[i];
  }

  /*! \brief Copy constructor. */
  template <class T> inline
  t3Matrix<T>::t3Matrix(const t3Matrix<T>& mat)
  {
    for (int i = 0; i < 9; i++) m[i] = mat.m[i];
  }

  /*! \brief Constructs a matrix from the given entries (row major). */
  template <class T> inline
  t3Matrix<T>::t3Matrix(T m0, T m3, T m6,
                        T m1, T m4, T m7,
                        T m2, T m5, T m8)
  {
    m[0] = m0;
    m[1] = m1;
    m[2] = m2;
    m[3] = m3;
    m[4] = m4;
    m[5] = m5;
    m[6] = m6;
    m[7] = m7;
    m[8] = m8;
  }

  /*! \brief Overwrites this matrix with the given entries (row major). */
  template <class T> inline
  void t3Matrix<T>::set(T m0, T m3, T m6,
                        T m1, T m4, T m7,
                        T m2, T m5, T m8)
  {
    m[0] = m0;
    m[1] = m1;
    m[2] = m2;
    m[3] = m3;
    m[4] = m4;
    m[5] = m5;
    m[6] = m6;
    m[7] = m7;
    m[8] = m8;
  }

  /*! \brief Scalar times matrix. */
  template <class T> inline
  t3Matrix<T> operator*(T f, const t3Matrix<T>& mat)
  {
    t3Matrix<T> out(mat);
    out *= f;
    return out;
  }

  /*! \brief Matrix times scalar. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::operator*(T f) const
  {
	t3Matrix<T> out(*this);
	out *= f;
	return out;
  }

  /*! \brief Matrix times scalar. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::operator/(T f) const
  {
	t3Matrix<T> out(*this);
	out /= f;
	return out;
  }

  /*! \brief Multiply this matrix with a scalar. */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::operator*=(const T& f)
  {
    for (int i = 0; i < 9; i++) m[i] *= f;
    return *this;
  }

  /*! \brief Divide this matrix by a scalar. */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::operator/=(const T& f)
  {
    for (int i = 0; i < 9; i++) m[i] /= f;
    return *this;
  }

  /*! \brief Matrix times vector. */
  template <class T> inline
  t3Vector<T> t3Matrix<T>::operator* (const t3Vector<T>& op) const
  {
    return t3Vector<T>(
             m[0]*op.x + m[3]*op.y + m[6]*op.z,
             m[1]*op.x + m[4]*op.y + m[7]*op.z,
             m[2]*op.x + m[5]*op.y + m[8]*op.z
           );
  }

  /*! \brief Matrix times vector. */
  template <class T> inline
  tRGB<T> t3Matrix<T>::operator* (const tRGB<T>& op) const
  {
    return tRGB<T>(
             m[0]*op.r + m[3]*op.g + m[6]*op.b,
             m[1]*op.r + m[4]*op.g + m[7]*op.b,
             m[2]*op.r + m[5]*op.g + m[8]*op.b
           );
  }

  /*! \brief Matrix times matrix. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::operator* (const t3Matrix<T>& op) const
  {

    return t3Matrix(m[0]*op.m[0] + m[3]*op.m[1] + m[6]*op.m[2],
                    m[0]*op.m[3] + m[3]*op.m[4] + m[6]*op.m[5],
                    m[0]*op.m[6] + m[3]*op.m[7] + m[6]*op.m[8],
                    m[1]*op.m[0] + m[4]*op.m[1] + m[7]*op.m[2],
                    m[1]*op.m[3] + m[4]*op.m[4] + m[7]*op.m[5],
                    m[1]*op.m[6] + m[4]*op.m[7] + m[7]*op.m[8],
                    m[2]*op.m[0] + m[5]*op.m[1] + m[8]*op.m[2],
                    m[2]*op.m[3] + m[5]*op.m[4] + m[8]*op.m[5],
                    m[2]*op.m[6] + m[5]*op.m[7] + m[8]*op.m[8]
                   );
  }

  /*! \brief Adds other matrix to this matrix. */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::operator+=(const t3Matrix<T>& op)
  {
    for (int i = 0; i < 9; i++) m[i] += op.m[i];
    return *this;
  }

  /*! \brief Subtracts other matrix from this matrix. */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::operator-=(const t3Matrix<T>& op)
  {
    *this += -op;
    return *this;
  }

  /*! \brief Matrix plus matrix. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::operator+(const t3Matrix<T>& op) const
  {
    t3Matrix<T> out(*this);
    return out += op;
  }

  /*! \brief Matrix minus matrix. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::operator-(const t3Matrix<T>& op) const
  {
    t3Matrix<T> out(*this);
    return out += -op;
  }

  /*! \brief Return additive inverse of this matrix. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::operator-() const
  {
    t3Matrix<T> out(*this);
    for (int i = 0; i < 9; i++) out[i] = -out[i];
    return out;
  }


  /*! \brief Right-multiply m to this matrix (*this = *this * m). */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::rmul(const t3Matrix<T>& m)
  {
    *this = *this * m;
    return *this;
  }

  /*! \brief Left-multiply m to this matrix (*this = m * *this). */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::lmul(const t3Matrix<T>& m)
  {
    *this = m * *this;
    return *this;
  }

  /*! \brief Return the trace of this matrix (\f$a_{11} + a_{22} + a_{33}\f$). */
  template <class T> inline
  T t3Matrix<T>::trace() const
  {
    return ( m[0] + m[4] + m[8] );
  }

  /*! \brief Return the determinant of this matrix. */
  template <class T> inline
  T t3Matrix<T>::det() const
  {
    return det3x3(m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]);
  }

  /*! \brief Return the adjugate of this matrix. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::adjugate() const
  {
    // transpose of cofactor matrix
    return t3Matrix(
             det2x2(m[4],m[7],m[5],m[8]), det2x2(m[5],m[8],m[3],m[6]), det2x2(m[3],m[6],m[4],m[7]),
             det2x2(m[7],m[1],m[8],m[2]), det2x2(m[8],m[2],m[6],m[0]), det2x2(m[6],m[0],m[7],m[1]),
             det2x2(m[1],m[4],m[2],m[5]), det2x2(m[2],m[5],m[0],m[3]), det2x2(m[0],m[3],m[1],m[4]));
  }

  /*! \brief Transpose this matrix.
   * Attention: Although innocent looking this is an inplace operation
   **/
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::transpose()
  {
    std::swap(m[1],m[3]);
    std::swap(m[2],m[6]);
    std::swap(m[5],m[7]);
    return *this;
  }

  /*! \brief Invert this matrix.
   * Attention: Although innocent looking this is an inplace operation
   **/
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::invert()
  {
    *this = (T(1)/det())*adjugate();
    return *this;
  }

  /*! \brief Overwrite this matrix with an identity matrix. */
  template <class T> inline
  t3Matrix<T>& t3Matrix<T>::loadIdentity ()
  {
    m[0] = T(1);
    m[1] = T(0);
    m[2] = T(0);
    m[3] = T(0);
    m[4] = T(1);
    m[5] = T(0);
    m[6] = T(0);
    m[7] = T(0);
    m[8] = T(1);
    return *this;
  }

  /*! \brief Retrieves the (not normalized) axis of rotation,
   *         assuming this matrix describes a rotation. */
  template <class T> inline
  t3Vector<T> t3Matrix<T>::getAxis() const
  {
    // gemaess Artin, "Algebra", Kapitel 4, Aufgabe 14
    float a0 = m[5] + m[7]; // (2,3) + (3,2)
    float a1 = m[2] + m[6]; // (1,3) + (3,1)
    float a2 = m[1] + m[3]; // (1,2) + (2,1)

    if (a0 == 0) return t3Vector<T>(T(1), T(0), T(0));
    else if (a1 == 0) return t3Vector<T>(T(0), T(1), T(0));
    else if (a2 == 0) return t3Vector<T>(T(0), T(0), T(1));
    else return t3Vector<T>(T(1)/a0, T(1)/a1, T(1)/a2);
  }

  /*! \brief Return the upper left 3x3 matrix from mat. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::extract(const t4Matrix<T>& mat)
  {

    return t3Matrix<T>(mat.m[0], mat.m[4], mat.m[8],
                       mat.m[1], mat.m[5], mat.m[9],
                       mat.m[2], mat.m[6], mat.m[10]);
  }

  /*! \brief Return a matrix representing a scaling by s. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::scale(const t3Vector<T>& s)
  {
    return t3Matrix<T>(
             s.x, T(0), T(0),
             T(0), s.y, T(0),
             T(0), T(0), s.z
           );
  }

  /*! Return a matrix that will rotate u into v. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::rotation(const t3Vector<T>& u, const t3Vector<T>& v)
  {
    T phi;
    T h;
    T lambda;
    t3Vector<T> w;

    w = u.cross(v);

    phi = u.dot(v);
    lambda = w.dot(w);

    if (lambda > 1e-10)
      h = ((T)1.0 - phi) / lambda;
    else
      h = lambda;

    T hxy = w.x * w.y * h;
    T hxz = w.x * w.z * h;
    T hyz = w.y * w.z * h;

    t3Matrix<T> out(phi + w.x * w.x * h, hxy + w.z, hxz - w.y,
                    hxy - w.z, phi + w.y * w.y * h, hyz + w.x,
                    hxz + w.y, hyz - w.x, phi + w.z * w.z * h);

    return out;
  }

  /*! \brief Return a matrix that rotates by specified angle (in degrees) around specified axis. */
  template <class T> inline
  t3Matrix<T> t3Matrix<T>::rotation(const t3Vector<T>& axis, float angle)
  {
    // formula copied form GL specification

    t3Vector<T> n(axis);
    n.normalize();
    // convert to radians
    angle *= (float)(3.1415927/180.);

    t3Matrix<T> s(0, -n.z, n.y,
                  n.z, 0, -n.x,
                  -n.y, n.x, 0);

    t3Matrix<T> nnt(n.x*n.x, n.x*n.y, n.x*n.z,
                    n.y*n.x, n.y*n.y, n.y*n.z,
                    n.z*n.x, n.z*n.y, n.z*n.z);

    return nnt + T(cos(angle))*(t3Matrix() - nnt) + T(sin(angle))*s;
  }

  template <class T> inline
  std::ostream& operator<< (std::ostream& os, const t3Matrix<T>& arg)
  {
    os << "[ " << std::setw(10) << arg[0] << "  " << std::setw(10) << arg[3] << "  " << std::setw(10) << arg[6] << " ]\n";
    os << "| " << std::setw(10) << arg[1] << "  " << std::setw(10) << arg[4] << "  " << std::setw(10) << arg[7] << " |\n";
    os << "| " << std::setw(10) << arg[2] << "  " << std::setw(10) << arg[5] << "  " << std::setw(10) << arg[8] << " |\n";
    return os;
  }

  template <class T> inline
  std::istream& operator>> ( std::istream& is, t3Matrix<T>& arg)
  {
    std::string dummy;
    is >> dummy >> arg[0] >> arg[3] >> arg[6] >> dummy;
    is >> dummy >> arg[1] >> arg[4] >> arg[7] >> dummy;
    is >> dummy >> arg[2] >> arg[5] >> arg[8] >> dummy;
    return is;
  }

  template <class T> inline
  t3Matrix<T> t3Matrix<T>::rotationX (T a)
  {
    return t3Matrix<T>(
             T(1), T(0),      T(0),
             T(0), T(cos(a)), T(-sin(a)),
             T(0), T(sin(a)), T(cos(a))
           );
  }

  template <class T> inline
  t3Matrix<T> t3Matrix<T>::rotationY (T a)
  {
    return t3Matrix<T>(
             T(cos(a)), T(0), T(-sin(a)),
             T(0),      T(1), T(0),
             T(sin(a)), T(0), T(cos(a))
           );
  }

  template <class T> inline
  t3Matrix<T> t3Matrix<T>::rotationZ (T a)
  {
    return t3Matrix<T>(
             T(cos(a)), T(-sin(a)), T(0),
             T(sin(a)), T(cos(a)),  T(0),
             T(0),      T(0),       T(1)
           );
  }

  typedef t3Matrix<float> f3Matrix;
  typedef t3Matrix<double> d3Matrix;


} // namespace gravis

#endif
