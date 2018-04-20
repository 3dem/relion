#ifndef __LIBGRAVIS_T_GRAY_A_H__
#define __LIBGRAVIS_T_GRAY_A_H__
/******************************************************************************
**        Title: tGray_A.h
**  Description: Represents an RGB+Alpha color tupel.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <iostream>

namespace gravis
{

  template <class T>
  class tGray_A
  {
    public:
      T g, a;

      typedef T scalar_type;

      tGray_A ()           : g(T(0)), a(T(1.0)) { }
      tGray_A (T _g)       : g(_g)  , a(T(1.0)) { }
      tGray_A (T _g, T _a) : g(_g)  , a(_a) { }

      void set (T _g)
      {
        g = _g;
      }

      void set (T _g, T _a)
      {
        g = _g;
        a = _a;
      }

      T grayValue () const
      {
        return g;
      }

      T minValue () const
      {
        return g;
      }

      T maxValue () const
      {
        return g;
      }

      /*! \brief All color components, including alpha are clamped to [0,1].
       *
       * \return self
       */
      tGray_A& clamp()
      {
        g = std::min(std::max(g, T(0)), T(1));
        return *this;
      }

      bool operator != (const tGray_A& c) const
      {
        return g != c.g || a != c.a;
      }

      bool operator == (const tGray_A& c) const
      {
        return g == c.g && a == c.a;
      }

      tGray_A& operator += (const tGray_A& c)
      {
        g += c.g;
        return *this;
      }

      tGray_A& operator += (const T gray)
      {
        g += gray;
        return *this;
      }

      tGray_A& operator -= (const tGray_A& c)
      {
        g -= c.g;
        return *this;
      }

      tGray_A& operator -= (const T gray)
      {
        g -= gray;
        return *this;
      }

      tGray_A& operator *= (const tGray_A& c)
      {
        g *= c.g;
        return *this;
      }

      tGray_A& operator *= (const float factor)
      {
        g *= factor;
        return *this;
      }

      tGray_A& operator /= (const tGray_A& c)
      {
        g /= c.g;
        return *this;
      }

      tGray_A& operator /= (const float factor)
      {
        g /= factor;
        return *this;
      }

      //! Unary minus
      inline
      tGray_A operator - () const
      {
        return tGray_A<T>(-g, a);
      };

      //! Addition of a scalar (analog to -=)
      inline
      tGray_A operator + (const T& c) const
      {
        return tGray_A<T>(g+c, a);
      };

      //! Subtraction of a scalar (analog to +=)
      inline
      tGray_A operator - (const T& c) const
      {
        return tGray_A<T>(g-c, a);
      };

      //! Multiplication of a scalar (analog to *=)
      inline
      tGray_A operator * (const T& c) const
      {
        return tGray_A<T>(g*c, a);
      };

      //! Division by a scalar (analog to /=)
      inline
      tGray_A operator / (const T& c) const
      {
        return tGray_A<T>(g/c, a);
      };
  };



  template <class T> inline
  tGray_A<T> operator+ (const tGray_A<T>& c1, const tGray_A<T>& c2)
  {
    tGray_A<T> result(c1);
    return (result += c2);
  }


  template <class T> inline
  tGray_A<T> operator- (const tGray_A<T>& c1, const tGray_A<T>& c2)
  {
    tGray_A<T> result(c1);
    return (result -= c2);
  }


  template <class T> inline
  tGray_A<T> operator* (const tGray_A<T>& c1, const tGray_A<T>& c2)
  {
    tGray_A<T> result(c1);
    return (result *= c2);
  }


  template <class T> inline
  tGray_A<T> operator* (const tGray_A<T>& c, T factor)
  {
    tGray_A<T> result(c);
    return (result *= factor);
  }


  template <class T> inline
  tGray_A<T> operator* (T factor, const tGray_A<T>& c)
  {
    tGray_A<T> result(c);
    return (result *= factor);
  }


  template <class T> inline
  tGray_A<T> operator / (const tGray_A<T>& c1, const tGray_A<T>& c2)
  {
    tGray_A<T> result(c1);
    return (result /= c2);
  }


  template <class T> inline
  tGray_A<T> operator / (const tGray_A<T>& c, T factor)
  {
    tGray_A<T> result(c);
    return (result /= factor);
  }


  template <class T> inline
  bool operator < (const tGray_A<T>& c1, const tGray_A<T>& c2)
  {
    return (c1.grayValue() < c2.grayValue());
  }


  template <class T> inline
  tGray_A<T> operator ! (const tGray_A<T>& c)
  {
    tGray_A<T> result = tGray_A<T>::White;
    return (result -= c);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tGray_A<T>& c)
  {
    os << "(" << c.g << " " << c.a << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tGray_A<unsigned char>& c)
  {
    os << "(" << (int)c.g << " " << (int)c.a << ")";
    return os;
  }

  typedef tGray_A<unsigned char>  bGray_A;
  typedef tGray_A<float>          fGray_A;
  typedef tGray_A<double>         dGray_A;

} /* Close Namespace "gravis" */

#endif
