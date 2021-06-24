#ifndef __LIBGRAVIS_T4VECTOR_H__
#define __LIBGRAVIS_T4VECTOR_H__
/******************************************************************************
 **        Title: t4Vector.h
 **  Description: Represents a four dimensional vector (3D+homogeneous comp.).
 **
 **       Author: Jean-Sebastien Pierrard, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace gravis
{

  template <class T> class t2Vector;
  template <class T> class t3Vector;

  template <class T>
  class t4Vector
  {
    public:
      T x, y, z, w;

      typedef T scalar_type;

      t4Vector () : x(T(0)), y(T(0)), z(T(0)), w(T(1)) { }
      explicit t4Vector (T _v) : x(_v), y(_v), z(_v), w(_v) { }
      t4Vector (T _x, T _y, T _z, T _w=T(1)) : x(_x), y(_y), z(_z), w(_w) { }
      /*! \brief Construct a 4D vector with w = 1. */
	  explicit t4Vector (const t3Vector<T>& vec) : x(vec.x), y(vec.y), z(vec.z), w(1.0) { }
	  explicit t4Vector (const t3Vector<T>& vec, T ww) : x(vec.x), y(vec.y), z(vec.z), w(ww) { }

      template <class T1>
      explicit t4Vector (const t4Vector<T1>& vec) : x(vec.x), y(vec.y), z(vec.z), w(vec.w) {}

      t4Vector (const t4Vector& vec)    : x(vec.x), y(vec.y), z(vec.z), w(vec.w) { }

      static t4Vector unitX ()
      {
        return t4Vector(T(1), T(0), T(0), T(1));
      }
      static t4Vector unitY ()
      {
        return t4Vector(T(0), T(1), T(0), T(1));
      }
      static t4Vector unitZ ()
      {
        return t4Vector(T(0), T(0), T(1), T(1));
      }

      void set (T _v)
      {
        x = y = z = _v;
        w = T(1);
      }

      void set (T _x, T _y, T _z, T _w=T(1))
      {
        x = _x;
        y = _y;
        z = _z;
        w = _w;
      }

      //! Beware: This is not the 2 norm but the square of the two norm.
      T norm2 () const
      {
        return (x*x + y*y + z*z + w*w);
      }

      //! \f$l_1\f$ Norm: \f$\sum_i |v_i|\f$
      T normL1 () const
      {
        return (std::abs(x) + std::abs(y) + std::abs(z) + std::abs(w));
      }

      //! \f$l_2\f$ Norm: \f$\sqrt{\sum_i |v_i|^2}\f$
      T normL2 () const
      {
        return sqrt(x*x + y*y + z*z + w*w);
      }

      //! \f$l_\infty\f$ Norm: \f$\max{ |v_i|\,|\, \forall i  }\f$
      T normLInf() const
      {
        return std::max(std::max(std::max(std::abs(x), std::abs(y)), std::abs(z)), std::abs(w));
      }

      void invert ()
      {
        x = -x;
        y = -y;
        z = -z;
        w = -w;
      }

      T dot (const t4Vector& arg) const
      {
        return (x*arg.x + y*arg.y + z*arg.z + w*arg.w);
      }

      void divideW ()
      {
        x /= w;
        y /= w;
        z /= w;
        w = T(1);
      }

      /*! \brief Return a 3D vector corresponding to this 4D vector.
       *
       * If the w coordinate is 0, the vector is considered a direction or displacement,
       * and (x,y,z) is returned. Otherwise, the vector is considered a point, and
       * (x/w, y/w, z/w) is returned.
       */
      t3Vector<T> toVector3() const
      {
        if (w == 0) return t3Vector<T>(x, y, z);
        else return t3Vector<T>(x/w, y/w, z/w);
      }
      
      t3Vector<T> xyz() const
	  {
		return t3Vector<T>(x, y, z);
	  }
	  
	  t2Vector<T> xy() const
	  {
		return t2Vector<T>(x, y);
	  }

      /*! \brief Return the euclidian norm of this 4D vector.
       *
       * Note, that there is no special treatment of the w-coordinate.
       * The result is simply \f$\sqrt{x^2+y^2+z^2+w^2}\f$.
       */
      T length () const
      {
        return T(::sqrt(x*x + y*y + z*z + w*w));
      }

      t4Vector& normalize (T f=T(1))
      {
        if (f == T(0)) set(T(0), T(0), T(0), T(0));
        T norm = length()/f;
        if (norm != T(0))
        {
          *this /= norm;
        }
        return *this;
      }

      const T& operator[] (int idx) const
      {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
        assert((idx >= 0) && (idx < 4));
#endif
        return (&x)[idx];
      }

      T& operator[] (int idx)
      {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
        assert((idx >= 0) && (idx < 4));
#endif
        return (&x)[idx];
      }

      bool operator == ( const t4Vector& arg ) const
      {
        return ( x == arg.x && y == arg.y && z == arg.z && w == arg.w);
      }

      bool operator != ( const t4Vector& arg ) const
      {
        return !(*this == arg);
      }

      t4Vector& operator += (const t4Vector& arg)
      {
        x += arg.x;
        y += arg.y;
        z += arg.z;
        w += arg.w;
        return *this;
      }

      t4Vector& operator -= (const t4Vector& arg)
      {
        x -= arg.x;
        y -= arg.y;
        z -= arg.z;
        w -= arg.w;
        return *this;
      }

      t4Vector& operator += (const T& scalar)
      {
        x += scalar;
        y += scalar;
        z += scalar;
        w += scalar;
        return *this;
      }

      t4Vector& operator -= (const T& scalar)
      {
        x -= scalar;
        y -= scalar;
        z -= scalar;
        w -= scalar;
        return *this;
      }

      t4Vector& operator *= (const T& arg)
      {
        x *= arg;
        y *= arg;
        z *= arg;
        w *= arg;
        return *this;
      }

      t4Vector& operator /= (const T& arg)
      {
        x /= arg;
        y /= arg;
        z /= arg;
        w /= arg;
        return *this;
      }

      //! Check if the entries of the other vector differ by less than epsilon.
      //  It is better to use this than to use operator== for comparision, if it is
      //  not the same vertex.
      bool isClose( const t4Vector& o, const T epsilon) const
      {
        return ((std::fabs(x-o.x) < epsilon) and (std::fabs(y-o.y) < epsilon) and (std::fabs(z-o.z) < epsilon) and (std::fabs(w-o.w) < epsilon));
      }

      static
      t4Vector normalize (const t4Vector& v1, T f=T(1))
      {
        return t4Vector(v1).normalize();
      }

      static
      T dot (const t4Vector& v1, const t4Vector& v2)
      {
        return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w);
      }
  };


  template <class T>
  inline
  t4Vector<T> operator + (const t4Vector<T>& v1, const t4Vector<T>& v2)
  {
    return t4Vector<T>(
             v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w
           );
  }


  template <class T>
  inline
  t4Vector<T> operator - (const t4Vector<T>& v1)
  {
    return t4Vector<T>(-v1.x, -v1.y, -v1.z, -v1.w);
  }


  template <class T>
  inline
  t4Vector<T> operator - (const t4Vector<T>& v1, const t4Vector<T>& v2)
  {
    return t4Vector<T>(
             v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w
           );
  }


  template <class T>
  inline
  t4Vector<T> operator + (const T& s, const t4Vector<T>& v2)
  {
    return t4Vector<T>(s + v2.x, s + v2.y, s + v2.z, s + v2.w);
  }

  template <class T>
  inline
  t4Vector<T> operator - (const T& s, const t4Vector<T>& v2)
  {
    return t4Vector<T>(s - v2.x, s - v2.y, s - v2.z, s - v2.w);
  }

  template <class T>
  inline
  t4Vector<T> operator + (const t4Vector<T>& v, const T& s)
  {
    return t4Vector<T>(v.x + s, v.y + s, v.z + s, v.w + s);
  }

  template <class T>
  inline
  t4Vector<T> operator - (const t4Vector<T>& v, const T& s)
  {
    return t4Vector<T>(v.x - s, v.y - s, v.z - s, v.w - s);
  }

  template <class T>
  inline
  t4Vector<T> operator * (T f, const t4Vector<T>& v)
  {
    return t4Vector<T>(f * v.x, f * v.y, f * v.z, f * v.w);
  }
  
  template <class T>
  inline
  t4Vector<T> operator * (const t4Vector<T>& v, T f)
  {
    return t4Vector<T>(f * v.x, f * v.y, f * v.z, f * v.w);
  }
  
  template <class T>
  inline
  t4Vector<T> operator / (const t4Vector<T>& v, T f)
  {
    return t4Vector<T>(v.x / f, v.y / f, v.z / f, v.w / f);
  }


  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const t4Vector<T>& arg)
  {
    os << "[" << arg.x << ", " << arg.y << ", " << arg.z << ", " << arg.w << "]";
    return os;
  }

  typedef t4Vector<float> f4Vector;
  typedef t4Vector<double> d4Vector;
  typedef t4Vector<int> i4Vector;

} /* Close Namespace "gravis" */

#endif
