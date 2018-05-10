#ifndef __LIBGRAVIS_T2VECTOR_H__
#define __LIBGRAVIS_T2VECTOR_H__
/******************************************************************************
**        Title: t2Vector.h
**  Description: Represents a two dimensional vector.
**
**       Author: Pascal Paysan, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace gravis
{

  template <class T>
  class t2Vector
  {
    public:
      T x, y;

      typedef T scalar_type;
      t2Vector() : x(T(0)), y(T(0)) { }
      explicit t2Vector(T _v) : x(_v), y(_v) { }
      t2Vector(T _x, T _y) : x(_x), y(_y) { }
      template <class T1>
      explicit t2Vector(const t2Vector<T1>& vec) : x(vec.x), y(vec.y) {}

      static t2Vector unitX()
      {
        return t2Vector(T(1), T(0));
      }
      static t2Vector unitY()
      {
        return t2Vector(T(0), T(1));
      }

      void set (T _v)
      {
        x = y= _v;
      }

      void set (T _x, T _y)
      {
        x = _x;
        y = _y;
      }

      T length () const
      {
        return T(::sqrt(x*x + y*y));
      }

      //! Beware: This is not the 2 norm but the square of the two norm.
      T norm2 () const
      {
        return (x*x + y*y);
      }

      //! \f$l_1\f$ Norm: \f$\sum_i |v_i|\f$
      T normL1 () const
      {
        return (std::abs(x) + std::abs(y));
      }

      //! \f$l_2\f$ Norm: \f$\sqrt{\sum_i |v_i|^2}\f$
      T normL2 () const
      {
        return sqrt(x*x + y*y);
      }

      //! \f$l_2\f$ Norm: \f$\sqrt{\sum_i |v_i|^2}\f$
      T normL2sqr () const
      {
        return x*x + y*y;
      }

      //! \f$l_\infty\f$ Norm: \f$\max{ |v_i|\,|\, \forall i  }\f$
      T normLInf() const
      {
        return std::max(std::abs(x), std::abs(y));
      }

      t2Vector& normalize (T f=1.0)
      {
        T norm = f / ::sqrt(x*x + y*y);
        x *= norm;
        y *= norm;
        return *this;
      }

      T dot (const t2Vector& arg) const
      {
        return (x*arg.x + y*arg.y);
      }

      const T& operator[] (int idx) const
      {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
        assert((idx >= 0) && (idx < 2));
#endif
        return (&x)[idx];
      }

      T& operator[] (int idx)
      {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
        assert((idx >= 0) && (idx < 2));
#endif
        return (&x)[idx];
      }

      bool operator == ( const t2Vector& arg ) const
      {
        return ( x == arg.x && y == arg.y);
      }

      bool operator != ( const t2Vector& arg ) const
      {
        return !(*this == arg);
      }

      t2Vector& operator += (const t2Vector& arg)
      {
        x += arg.x;
        y += arg.y;
        return *this;
      }

      t2Vector& operator -= (const t2Vector& arg)
      {
        x -= arg.x;
        y -= arg.y;
        return *this;
      }

      t2Vector& operator += (const T& scalar)
      {
        x += scalar;
        y += scalar;
        return *this;
      }

      t2Vector& operator -= (const T& scalar)
      {
        x -= scalar;
        y -= scalar;
        return *this;
      }

      t2Vector& operator *= (const T& arg)
      {
        x *= arg;
        y *= arg;
        return *this;
      }

      t2Vector& operator /= (const T& arg)
      {
        x /= arg;
        y /= arg;
        return *this;
      }

      //! Check if the entries of the other vector differ by less than epsilon.
      //  It is better to use this than to use operator== for comparision, if it is
      //  not the same vertex.
      bool isClose( const t2Vector& o, const T epsilon) const
      {
        return ((std::fabs(x-o.x) < epsilon) and (std::fabs(y-o.y) < epsilon));
      }

      static
      t2Vector normalize (const t2Vector& v1, T f=1.0f)
      {
        T norm = f / T(::sqrt(v1.x*v1.x + v1.y*v1.y));
        return t2Vector(v1.x * norm, v1.y * norm);
      }

      t2Vector operator / (const T& arg) const
      {
        t2Vector r(*this);
        r /= arg;
        return r;
      }

  };


  template <class T>
  inline
  t2Vector<T> operator + (const t2Vector<T>& v1, const t2Vector<T>& v2)
  {
    return t2Vector<T>(v1.x + v2.x, v1.y + v2.y);
  }


  template <class T>
  inline
  t2Vector<T> operator - (const t2Vector<T>& v1)
  {
    return t2Vector<T>(-v1.x, -v1.y);
  }


  template <class T>
  inline
  t2Vector<T> operator - (const t2Vector<T>& v1, const t2Vector<T>& v2)
  {
    return t2Vector<T>(v1.x - v2.x, v1.y - v2.y);
  }

  template <class T>
  inline
  t2Vector<T> operator + (const T& s, const t2Vector<T>& v2)
  {
    return t2Vector<T>(s + v2.x, s + v2.y);
  }

  template <class T>
  inline
  t2Vector<T> operator - (const T& s, const t2Vector<T>& v2)
  {
    return t2Vector<T>(s - v2.x, s - v2.y);
  }

  template <class T>
  inline
  t2Vector<T> operator + (const t2Vector<T>& v, const T& s)
  {
    return t2Vector<T>(v.x + s, v.y + s);
  }

  template <class T>
  inline
  t2Vector<T> operator - (const t2Vector<T>& v, const T& s)
  {
    return t2Vector<T>(v.x - s, v.y - s);
  }

  template <class T>
  inline
  t2Vector<T> operator * (T f, const t2Vector<T>& v)
  {
    return t2Vector<T>(f * v.x, f * v.y);
  }

  template <class T>
  inline
  t2Vector<T> operator * (const t2Vector<T>& v, const T& f)
  {
    return t2Vector<T>(f * v.x, f * v.y);
  }

  template <class T>
  inline
  t2Vector<T> operator * (const t2Vector<T>& v, const t2Vector<T>& f)
  {
    return t2Vector<T>(v.x * f.x, v.y * f.y);
  }

  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const t2Vector<T>& arg)
  {
    os << "[" << arg.x << ", " << arg.y << "]";
    return os;
  }

  template <class T>
  inline
  T dot (const t2Vector<T>& v1, const t2Vector<T>& v2)
  {
    return (v1.x*v2.x + v1.y*v2.y);
  }

  // Inverse of operator<<
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, t2Vector<T>& arg)
  {
    char c = ' ';
    is >> c;
    if (is.eof())
      return is;
    if (c != '[')
      throw std::runtime_error("Vector should start with an opening [");
    std::stringstream values;
    int v = 0;
    while ((is >> c) && (c != ']'))
    {
      if (c == ',')
      {
        v++;
        if (v >= 2)
          throw std::runtime_error("Vector contains more than three elements");
        values << " ";
      }
      else if (c != ' ')
        values << c;
    }
    if (c != ']')
    {
      throw std::runtime_error("Vector should end with a ]");
    }
    values >> arg.x >> arg.y;
    return is;
  }

  typedef t2Vector<float > f2Vector;
  typedef t2Vector<double> d2Vector;
  typedef t2Vector<int> i2Vector;

} /* Close Namespace "gravis" */

#endif
