#ifndef __LIBGRAVIS_T_RGB_A_H__
#define __LIBGRAVIS_T_RGB_A_H__
/******************************************************************************
**        Title: tRGB_A.h
**  Description: Represents an RGB+Alpha color tupel.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include "tRGB.h"
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace gravis
{

  template <class T>
  class tRGB_A
  {
    public:
      T r, g, b, a;

      tRGB_A () : r(T(0)), g(T(0)), b(T(0)), a(T(1.0)) { }
      tRGB_A (T _r, T _g, T _b, T _a=T(1.0)) : r(_r), g(_g), b(_b), a(_a) { }
      tRGB_A (T gray) : r(gray), g(gray), b(gray), a(T(1.0)) { }
      tRGB_A (T gray, T alpha) : r(gray), g(gray), b(gray), a(alpha) { }
      explicit tRGB_A (const tRGBA<T>& c) : r(c.r), g(c.g), b(c.b), a(c.a) { }
      explicit tRGB_A (const tRGB<T>& c, T _a=T(1.0)) : r(c.r), g(c.g), b(c.b), a(_a) { }

      void set (T _r, T _g, T _b, T _a)
      {
        r = _r;
        g = _g;
        b = _b;
        a = _a;
      }

      void set (T _r, T _g, T _b)
      {
        r = _r;
        g = _g;
        b = _b;
      }

      void set (T gray)
      {
        r = gray;
        b = gray;
        g = gray;
      }

      void set (T gray, T alpha)
      {
        r = gray;
        b = gray;
        g = gray;
        a = alpha;
      }

      T grayValue () const
      {
        return (T)(0.30*r + 0.59*g + 0.11*b);
      }

      T minValue () const
      {
        if (r < g)
        {
          if (r < b) return r;
          else return b;
        }
        else
        {
          if (g < b) return g;
          else return b;
        }
      }

      T maxValue () const
      {
        if (r > g)
        {
          if (r > b) return r;
          else return b;
        }
        else
        {
          if (g > b) return g;
          else return b;
        }
      }

      /*! \brief All color components, including alpha are clamped to [0,1].
       *
       * \return self
       */
      tRGB_A& clamp()
      {
        r = std::min(std::max(r, T(0)), T(1));
        g = std::min(std::max(g, T(0)), T(1));
        b = std::min(std::max(b, T(0)), T(1));
        return *this;
      }

      bool operator != (const tRGB_A& c) const
      {
        return r != c.r || g != c.g || b != c.b || a != c.a;
      }

      bool operator == (const tRGB_A& c) const
      {
        return r == c.r && g == c.g && b == c.b && a == c.a;
      }

      tRGB_A& operator += (const tRGB_A& c)
      {
        r += c.r;
        g += c.g;
        b += c.b;
        return *this;
      }

      tRGB_A& operator += (const T gray)
      {
        r += gray;
        g += gray;
        b += gray;
        return *this;
      }

      tRGB_A& operator -= (const tRGB_A& c)
      {
        r -= c.r;
        g -= c.g;
        b -= c.b;
        return *this;
      }

      tRGB_A& operator -= (const T gray)
      {
        r -= gray;
        g -= gray;
        b -= gray;
        return *this;
      }

      tRGB_A& operator *= (const tRGB_A& c)
      {
        r *= c.r;
        g *= c.g;
        b *= c.b;
        return *this;
      }

      tRGB_A& operator *= (const float factor)
      {
        r *= factor;
        g *= factor;
        b *= factor;
        return *this;
      }

      tRGB_A& operator /= (const tRGB_A& c)
      {
        r /= c.r;
        g /= c.g;
        b /= c.b;
        return *this;
      }

      tRGB_A& operator /= (const float factor)
      {
        r /= factor;
        g /= factor;
        b /= factor;
        return *this;
      }

      //! Unary minus
      inline
      tRGB_A operator - () const
      {
        return tRGB_A<T>(-r, -g, -b, a);
      };

      //! Addition of a scalar (analog to -=)
      inline
      tRGB_A operator + (const T& c) const
      {
        return tRGB_A<T>(r+c, g+c, b+c, a);
      };

      //! Subtraction of a scalar (analog to +=)
      inline
      tRGB_A operator - (const T& c) const
      {
        return tRGB_A<T>(r-c, g-c, b-c, a);
      };

      //! Multiplication of a scalar (analog to *=)
      inline
      tRGB_A operator * (const T& c) const
      {
        return tRGB_A<T>(r*c, g*c, b*c, a);
      };

      //! Division by a scalar (analog to /=)
      inline
      tRGB_A operator / (const T& c) const
      {
        return tRGB_A<T>(r/c, g/c, b/c, a);
      };
  };



  template <class T> inline
  tRGB_A<T> operator+ (const tRGB_A<T>& c1, const tRGB_A<T>& c2)
  {
    tRGB_A<T> result(c1);
    return (result += c2);
  }


  template <class T> inline
  tRGB_A<T> operator- (const tRGB_A<T>& c1, const tRGB_A<T>& c2)
  {
    tRGB_A<T> result(c1);
    return (result -= c2);
  }


  template <class T> inline
  tRGB_A<T> operator* (const tRGB_A<T>& c1, const tRGB_A<T>& c2)
  {
    // tRGB_A<T> result(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b, c1.a * c2.a);
    tRGB_A<T> result(c1);
    result *= c2;
    return result;
  }


  template <class T> inline
  tRGB_A<T> operator* (const tRGB_A<T>& c, T factor)
  {
    // 	tRGB_A<T> result(c.r * factor, c.g * factor, c.b * factor, c.a);
    tRGB_A<T> result(c);
    return (result *= factor);
  }


  template <class T> inline
  tRGB_A<T> operator* (T factor, const tRGB_A<T>& c)
  {
    // 	tRGB_A<T> result(c.r * factor, c.g * factor, c.b * factor, c.a);
    tRGB_A<T> result(c);
    return (result *= factor);
  }


  template <class T> inline
  tRGB_A<T> operator / (const tRGB_A<T>& c1, const tRGB_A<T>& c2)
  {
    tRGB_A<T> result(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b);
    return result;
  }


  template <class T> inline
  tRGB_A<T> operator / (const tRGB_A<T>& c, T factor)
  {
    tRGB_A<T> result(c.r / factor, c.g / factor, c.b / factor);
    return result;
  }


  template <class T> inline
  bool operator < (const tRGB_A<T>& c1, const tRGB_A<T>& c2)
  {
    T gray1 = c1.grayValue();
    T gray2 = c2.grayValue();
    return (gray1 < gray2);
  }


  template <class T> inline
  tRGB_A<T> operator ! (const tRGB_A<T>& c)
  {
    tRGB_A<T> result = tRGB_A<T>::White;
    return (result -= c);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tRGB_A<T>& c)
  {
    os << "(" << c.r << " " << c.g << " " << c.b << " " << c.a << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tRGB_A<unsigned char>& c)
  {
    os << "(" << (int)c.r << " " << (int)c.g
       << " " << (int)c.b << " " << (int)c.a << ")";
    return os;
  }
  
    // Inverse of operator<<
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, tRGB_A<T>& arg)
  {
    char c = ' ';
    is >> c;
    if (is.eof())
      return is;
    if (c != '(')
      throw std::runtime_error("tRGB_A should start with an opening (");
    std::stringstream values;
    int v = 0;
    while ((is >> c) && (c != ')'))
    {
      if (c == ' ')
      {
        v++;
        if (v >= 4)
          throw std::runtime_error("tRGB_A contains more than four elements");
        values << " ";
      }
      else if (c != ' ')
        values << c;
    }
    if (c != ')')
    {
      throw std::runtime_error("tRGB_A should end with a )");
    }
    if ( v < 3 )
    {
      throw std::runtime_error("tRGB_A has not enough color values");
    }
    values >> arg.r >> arg.g >> arg.b >> arg.a;
    return is;
  }

  typedef tRGB_A<unsigned char>  bRGB_A;
  typedef tRGB_A<float>          fRGB_A;
  typedef tRGB_A<double>         dRGB_A;

} /* Close Namespace "gravis" */

#endif
