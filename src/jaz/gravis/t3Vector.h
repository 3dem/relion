#ifndef __LIBGRAVIS_T3VECTOR_H__
#define __LIBGRAVIS_T3VECTOR_H__
/******************************************************************************
**        Title: t3Vector.h
**  Description:
**
******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace gravis
{

  template <class T> class t4Vector;
  template <class T> class t2Vector;
  template <class T> class tRGBA;

  template <class T>
  class t3Vector
  {
    public:

      typedef T scalar_type;

      t3Vector () : x(T(0)), y(T(0)), z(T(0)) {}
      explicit t3Vector (T _v) : x(_v), y(_v), z(_v) {}
      t3Vector (T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
      t3Vector (const t3Vector& vec) : x(vec.x), y(vec.y), z(vec.z) {}
      // Generalized Copy Constructor (allows e.g. conversion from double to float)
      template <class T1>
      explicit
      t3Vector (const t3Vector<T1>& vec) : x(vec.x), y(vec.y), z(vec.z) {}
      // Initialization from Array
      t3Vector (const T* vecd) : x(vecd[0]), y(vecd[1]), z(vecd[2]) {}

      explicit
      t3Vector (const t4Vector<T>& vec) : x(vec.x/vec.w), y(vec.y/vec.w), z(vec.z/vec.w) {}     //rk

      /*
        Deprecated. Should use t4Vector::toVector3() which has more appropriate logic
        - this here makes only sense under the (wrong) assumption that all 4D vectors represent finite points
        - but homogeneous coordinates can also be directions, displacements or points at infinity,
        all with w == 0. - in any case, too much logic for a cast - BAD (mk)

      */

      static t3Vector unitX ()
      {
        return t3Vector(T(1), T(0), T(0));
      }
      static t3Vector unitY ()
      {
        return t3Vector(T(0), T(1), T(0));
      }
      static t3Vector unitZ ()
      {
        return t3Vector(T(0), T(0), T(1));
      }


      void set (T _v)
      {
        x = y = z = _v;
      }

      void set (T _x, T _y, T _z)
      {
        x = _x;
        y = _y;
        z = _z;
      }

      T length () const
      {
        return ::sqrt(x*x + y*y + z*z);
      }

      //! Beware: This is not the 2 norm but the square of the two norm.
      T norm2 () const
      {
        return (x*x + y*y + z*z);
      }

      //! Squared L2 Norm
      T normL2sqr () const
      {
        return (x*x + y*y + z*z);
      }

      //! \f$l_1\f$ Norm: \f$\sum_i |v_i|\f$
      T normL1 () const
      {
        return (std::fabs(x) + std::fabs(y) + std::fabs(z));
      }

      //! \f$l_2\f$ Norm: \f$\sqrt{\sum_i |v_i|^2}\f$
      T normL2 () const
      {
        return sqrt(x*x + y*y + z*z);
      }

      //! \f$l_\infty\f$ Norm: \f$\max{ |v_i|\,|\, \forall i  }\f$
      T normLInf() const
      {
        return std::max(std::max(std::abs(x), std::abs(y)), std::abs(z));
      }

      T sum() const
      {
        return x + y + z;
      }

      t3Vector<T> findOrthogonal() const;

      void invert ()
      {
        x = -x;
        y = -y;
        z = -z;
      }

      T dot (const t3Vector& arg) const
      {
        return (x*arg.x + y*arg.y + z*arg.z);
      }

      t3Vector cross (const t3Vector& arg) const
      {
        return t3Vector<T>(
                 y*arg.z - z*arg.y,
                 z*arg.x - x*arg.z,
                 x*arg.y - y*arg.x
               );
      }

      /**
       * Inplace normalization
       **/
      t3Vector& normalize (T f=T(1))
      {
        if (f == T(0)) set(T(0), T(0), T(0));
        T norm = length()/f;
        if (norm != T(0))
        {
          *this /= norm;
        }
        return *this;
      }

      /*! \brief Component wise multiplication (matlab ".*"). */
      t3Vector cmul(const t3Vector& v) const
      {
        return t3Vector(x * v.x, y * v.y, z * v.z);
      }

      /*! \brief Component wise division (matlab "./"). */
      t3Vector cdiv(const t3Vector& v) const
      {
        return t3Vector(x / v.x, y / v.y, z / v.z);
      }

      /*! \brief Interpolate three values, using this vector as barycentric coordinates. */
      template <class Value>
      Value interpolate(const Value& a, const Value& b, const Value& c) const
      {
        return x * a + y * b + z * c;
      }
      
      t2Vector<T> xy() const
	  {
		return t2Vector<T>(x, y);
	  }

      const T& operator[] (int idx) const
      {
        return (&x)[idx];
      }

      T& operator[] (int idx)
      {
        return (&x)[idx];
      }

      bool operator == ( const t3Vector& arg ) const
      {
        return ( x == arg.x && y == arg.y && z == arg.z );
      }

      bool operator != ( const t3Vector& arg ) const
      {
        return !(*this == arg);
      }

      t3Vector& operator += (const t3Vector& arg)
      {
        x += arg.x;
        y += arg.y;
        z += arg.z;
        return *this;
      }

      t3Vector& operator -= (const t3Vector& arg)
      {
        x -= arg.x;
        y -= arg.y;
        z -= arg.z;
        return *this;
      }

      t3Vector& operator += (const T& scalar)
      {
        x += scalar;
        y += scalar;
        z += scalar;
        return *this;
      }

      t3Vector& operator -= (const T& scalar)
      {
        x -= scalar;
        y -= scalar;
        z -= scalar;
        return *this;
      }

      t3Vector& operator *= (T arg)
      {
        x *= arg;
        y *= arg;
        z *= arg;
        return *this;
      }

      t3Vector operator * (T arg) const
      {
        return t3Vector(x * arg, y * arg, z * arg);
      }

      t3Vector& operator /= (T arg)
      {
        x /= arg;
        y /= arg;
        z /= arg;
        return *this;
      }

      t3Vector operator / (T arg) const
      {
        return t3Vector(x / arg, y / arg, z / arg);
      }

      T dist2( const t3Vector& v ) const
      {
        return ((x-v.x)*(x-v.x)+(y-v.y)*(y-v.y)+(z-v.z)*(z-v.z));
      }

      T dist( const t3Vector& v ) const
      {
        return ::sqrt( dist2( v ) );
      }

      //! Check if the entries of the other vector differ by less than epsilon.
      //  It is better to use this than to use operator== for comparision, if it is
      //  not the same vertex.
      bool isClose( const t3Vector& o, const T epsilon) const
      {
        return ((std::fabs(x-o.x) < epsilon) and (std::fabs(y-o.y) < epsilon) and (std::fabs(z-o.z) < epsilon));
      }

      static
      t3Vector normalize (const t3Vector& v1, T f=T(1))
      {
        return t3Vector(v1).normalize(f);
      }

      static
      T dot (const t3Vector& v1, const t3Vector& v2)
      {
        return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
      }

      static
      t3Vector cross (const t3Vector& v1, const t3Vector& v2)
      {
        return t3Vector(
                 v1.y*v2.z - v1.z*v2.y,
                 v1.z*v2.x - v1.x*v2.z,
                 v1.x*v2.y - v1.y*v2.x
               );
      }

    public:
      T x, y, z;
  };

  /*! \brief Returns minimal components of two vectors.
   *
   * This is useful for quick and dirty calculations of boundings boxes
   * for a set of vectors.
   *
   * TODO: This should be a static function of t3Vector
   */
  template <class T> inline
  t3Vector<T> lowerBound(const t3Vector<T> v1, const t3Vector<T> v2)
  {
    return t3Vector<T>(min(v1.x, v2.x), min(v1.y, v2.y), min(v1.z, v2. z));
  }

  /*! \brief Returns maximal components of two vectors.
   *
   * This is useful for quick and dirty calculations of boundings boxes
   * for a set of vectors.
   *
   * TODO: This should be a static function of t3Vector
   */
  template <class T> inline
  t3Vector<T> upperBound(const t3Vector<T> v1, const t3Vector<T> v2)
  {
    return t3Vector<T>(max(v1.x, v2.x), max(v1.y, v2.y), max(v1.z, v2. z));
  }

  /*! \brief Returns a vector orthogonal to this vector.
   *
   * E.g. for ||y|| < ||x|| and ||y|| < ||z||, the returned
   * vector is (z, 0, -x).
   */
  template <class T> inline
  t3Vector<T> t3Vector<T>::findOrthogonal() const
  {
    if (std::abs(y) < std::abs(z))
    {
      // y < z
      if (std::abs(x) < std::abs(y))
      {
        // x smallest
        return t3Vector<T>(0, z, -y);
      }
      else
      {
        // y smallest
        return t3Vector<T>(z, 0, -x);
      }
    }
    else
    {
      // z < y
      if (std::abs(x) < std::abs(z))
      {
        // x smallest
        return t3Vector<T>(0, z, -y);
      }
      else
      {
        // z smallest
        return t3Vector<T>(y, -x, 0);
      }
    }
  }


  template <class T> inline
  t3Vector<T> operator ~ (const t3Vector<T>& v1)
  {
    return t3Vector<T>(-v1.x, -v1.y, -v1.z);
  }


  template <class T>
  inline
  t3Vector<T> operator + (const t3Vector<T>& v1, const t3Vector<T>& v2)
  {
    return t3Vector<T>(
             v1.x + v2.x, v1.y + v2.y, v1.z + v2.z
           );
  }


  template <class T>
  inline
  t3Vector<T> operator - (const t3Vector<T>& v1)
  {
    return t3Vector<T>(-v1.x, -v1.y, -v1.z);
  }


  template <class T>
  inline
  t3Vector<T> operator + (const T& s, const t3Vector<T>& v2)
  {
    return t3Vector<T>(s + v2.x, s + v2.y, s + v2.z);
  }

  template <class T>
  inline
  t3Vector<T> operator - (const T& s, const t3Vector<T>& v2)
  {
    return t3Vector<T>(s - v2.x, s - v2.y, s - v2.z);
  }

  template <class T>
  inline
  t3Vector<T> operator + (const t3Vector<T>& v, const T& s)
  {
    return t3Vector<T>(v.x + s, v.y + s, v.z + s);
  }

  template <class T>
  inline
  t3Vector<T> operator - (const t3Vector<T>& v, const T& s)
  {
    return t3Vector<T>(v.x - s, v.y - s, v.z - s);
  }

  template <class T>
  inline
  t3Vector<T> operator - (const t3Vector<T>& v1, const t3Vector<T>& v2)
  {
    return t3Vector<T>(
             v1.x - v2.x, v1.y - v2.y, v1.z - v2.z
           );
  }

  template <class T>
  inline
  t3Vector<T> operator * (T f, const t3Vector<T>& v)
  {
    return t3Vector<T>(f * v.x, f * v.y, f * v.z);
  }

  template <class T>
  inline
  t3Vector<T> operator * (const t3Vector<T>& v, T f)
  {
    return t3Vector<T>(f * v.x, f * v.y, f * v.z);
  }

  template <class T>
  inline
  t3Vector<T> operator / (const t3Vector<T>& v, T f)
  {
    return t3Vector<T>( v.x/f, v.y/f, v.z/f );
  }

  template <class T>
  inline
  bool operator < (const t3Vector<T>& v1, const t3Vector<T>& v2)
  {
    return ((v1.x  < v2.x) || ((v1.x == v2.x) && (v1.y  < v2.y)) || ((v1.x == v2.x) && (v1.y == v2.y) && (v1.z < v2.z)));
  }

  template <class T>
  inline
  bool operator > (const t3Vector<T>& v1, const t3Vector<T>& v2)
  {
    return (!((v1==v2) || (v1<v2)));
  }

  template <class T>
  inline
  T dot (const t3Vector<T>& v1, const t3Vector<T>& v2)
  {
    return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
  }

  template <class T>
  inline
  t3Vector<T> cross (const t3Vector<T>& v1, const t3Vector<T>& v2)
  {
    return t3Vector<T>(
             v1.y*v2.z - v1.z*v2.y,
             v1.z*v2.x - v1.x*v2.z,
             v1.x*v2.y - v1.y*v2.x
           );
  }

  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const t3Vector<T>& arg)
  {
	os << "[" << arg.x << ", " << arg.y << ", " << arg.z << "]";
    return os;
  }

  // Inverse of operator<<
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, t3Vector<T>& arg)
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
        if (v >= 3)
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
    values >> arg.x >> arg.y >> arg.z;
    return is;
  }

  typedef t3Vector<float> f3Vector;
  typedef t3Vector<double> d3Vector;
  typedef t3Vector<int> i3Vector;

} // namespace gravis

#endif
