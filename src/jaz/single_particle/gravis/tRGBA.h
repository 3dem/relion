#ifndef __LIBGRAVIS_T_RGBA_H__
#define __LIBGRAVIS_T_RGBA_H__
/******************************************************************************
**        Title: tRGBA.h
**  Description: Represents an RGB+Alpha color tupel.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Brian Amberg, 2005-2006
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
  class tRGBA
  {
    public:
      T r, g, b, a;

      typedef T scalar_type;
      //! Default constructs a black, translucent pixel
      tRGBA () : r(T(0)), g(T(0)), b(T(0)), a(T(0)) { }
      tRGBA (const T& r, const T& g, const T& b, const T& a=T(1)) : r(r), g(g), b(b), a(a) { }
      tRGBA (const T& gray, const T& alpha=T(1)) : r(gray), g(gray), b(gray), a(alpha) { }
      explicit tRGBA (const tRGB<T>& c, const T& a=T(1)) : r(c.r), g(c.g), b(c.b), a(a) { }

      void set (T _r, T _g, T _b, T _a)
      {
        r = _r;
        g = _g;
        b = _b;
        a = _a;
      }

      void set (T gray)
      {
        r = gray;
        b = gray;
        g = gray;
        a = T(1);
      }

      /*!
       * Conversion to a gray pixel
       *
       * TODO: This should be put in an external conversion file, together with cie, hsv, etc...
       **/
      T grayValue () const
      {
        return T(0.30*r + 0.59*g + 0.11*b);
      }

      T intensity () const
      {
        return grayValue();
      }

      /*!
       * Return minimum of the tupel, ignoring the alpha channel.
       *
       * TODO: Is this really necessary in here. It could be a utility function.
       **/
      T const& minValue () const
      {
        return std::min(std::min(r, g), b);
      }

      /*!
       * Return maximum of the tupel, ignoring the alpha channel.
       *
       * TODO: Is this really necessary in here. It could be a utility function.
       **/
      T const& maxValue () const
      {
        return std::max(std::max(r, g), b);
      }

      const T& operator [] (const size_t& i) const
      {
        return (&r)[i];
      }
      T& operator [] (const size_t& i)
      {
        return (&r)[i];
      }

      /*!
       * \brief All color components, including alpha are clamped to [0,1]. This function works inplace.
       *
       * \return self
       */
      tRGBA& clamp()
      {
        r = std::min(std::max(r, T(0)), T(1));
        g = std::min(std::max(g, T(0)), T(1));
        b = std::min(std::max(b, T(0)), T(1));
        a = std::min(std::max(a, T(0)), T(1));
        return *this;
      }

      bool operator != (const tRGBA& c) const
      {
        return r != c.r || g != c.g || b != c.b || a != c.a;
      }

      bool operator == (const tRGBA& c) const
      {
        return r == c.r && g == c.g && b == c.b && a == c.a;
      }

      tRGBA& operator += (const tRGBA& c)
      {
        r += c.r;
        g += c.g;
        b += c.b;
        a += c.a;
        return *this;
      }

//      tRGBA& operator += (const T gray)
//      {
//        r += gray;
//        g += gray;
//        b += gray;
//        return *this;
//      }

      tRGBA& operator -= (const tRGBA& c)
      {
        r -= c.r;
        g -= c.g;
        b -= c.b;
        a -= c.a;
        return *this;
      }

//      tRGBA& operator -= (const T gray)
//      {
//        r -= gray;
//        g -= gray;
//        b -= gray;
//        return *this;
//      }

      tRGBA& operator *= (const tRGBA& c)
      {
        r *= c.r;
        g *= c.g;
        b *= c.b;
        a *= c.a;
        return *this;
      }

      tRGBA& operator *= (const float factor)
      {
        r *= factor;
        g *= factor;
        b *= factor;
        a *= factor;
        return *this;
      }

      tRGBA& operator /= (const tRGBA& c)
      {
        r /= c.r;
        g /= c.g;
        b /= c.b;
        a /= c.a;
        return *this;
      }

      tRGBA& operator /= (const float factor)
      {
        r /= factor;
        g /= factor;
        b /= factor;
        a /= factor;
        return *this;
      }

      //! Unary minus
      inline
      tRGBA operator - () const
      {
        return tRGBA<T>(-r, -g, -b, -a);
      };

      //! Addition of a scalar (analog to -=)
//      inline
//      tRGBA operator + (const T& c) const
//      {
//        return tRGBA<T>(r+c, g+c, b+c, a);
//      };

      //! Subtraction of a scalar (analog to +=)
//      inline
//      tRGBA operator - (const T& c) const
//      {
//        return tRGBA<T>(r-c, g-c, b-c, a);
//      };

      //! Multiplication of a scalar (analog to *=)
      inline
      tRGBA operator * (const T& c) const
      {
        return tRGBA<T>(r*c, g*c, b*c, a*c);
      };

      //! Division by a scalar (analog to /=)
      inline
      tRGBA operator / (const T& c) const
      {
        return tRGBA<T>(r/c, g/c, b/c, a/c);
      };
  };



  template <class T> inline
  tRGBA<T> operator+ (const tRGBA<T>& c1, const tRGBA<T>& c2)
  {
    tRGBA<T> result(c1);
    return (result += c2);
  }


  template <class T> inline
  tRGBA<T> operator- (const tRGBA<T>& c1, const tRGBA<T>& c2)
  {
    tRGBA<T> result(c1);
    return (result -= c2);
  }


  template <class T> inline
  tRGBA<T> operator* (const tRGBA<T>& c1, const tRGBA<T>& c2)
  {
    tRGBA<T> result(c1);
    result *= c2;
    return result;
  }


  template <class T> inline
  tRGBA<T> operator* (const tRGBA<T>& c, T factor)
  {
    tRGBA<T> result(c);
    return (result *= factor);
  }


  template <class T> inline
  tRGBA<T> operator* (T factor, const tRGBA<T>& c)
  {
    tRGBA<T> result(c);
    return (result *= factor);
  }


  template <class T> inline
  tRGBA<T> operator / (const tRGBA<T>& c1, const tRGBA<T>& c2)
  {
    tRGBA<T> result(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b);
    return result;
  }


  template <class T> inline
  tRGBA<T> operator / (const tRGBA<T>& c, T factor)
  {
    tRGBA<T> result(c.r / factor, c.g / factor, c.b / factor);
    return result;
  }


  template <class T> inline
  bool operator < (const tRGBA<T>& c1, const tRGBA<T>& c2)
  {
    T gray1 = c1.grayValue();
    T gray2 = c2.grayValue();
    return (gray1 < gray2);
  }


  template <class T> inline
  tRGBA<T> operator ! (const tRGBA<T>& c)
  {
    tRGBA<T> result = tRGBA<T>::White;
    return (result -= c);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tRGBA<T>& c)
  {
    os << "(" << c.r << " " << c.g << " " << c.b << " " << c.a << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tRGBA<unsigned char>& c)
  {
    os << "(" << (int)c.r << " " << (int)c.g
       << " " << (int)c.b << " " << (int)c.a << ")";
    return os;
  }
  
  // Inverse of operator<<
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, tRGBA<T>& arg)
  {
    char c = ' ';
    is >> c;
    if (is.eof())
      return is;
    if (c != '(')
      throw std::runtime_error("tRGBA should start with an opening (");
    std::stringstream values;
    int v = 0;
    while ((is >> c) && (c != ')'))
    {
      if (c == ',')
      {
        v++;
        if (v >= 4)
          throw std::runtime_error("tRGBA contains more than four elements");
        values << " ";
      }
      else if (c != ' ')
        values << c;
    }
    if (c != ')')
    {
      throw std::runtime_error("tRGBA should end with a )");
    }
    if ( v < 3 )
    {
      throw std::runtime_error("tRGBA has not enough color values");
    }
    values >> arg.r >> arg.g >> arg.b >> arg.a;
    return is;
  }


  typedef tRGBA<unsigned char>  bRGBA;
  typedef tRGBA<float>          fRGBA;
  typedef tRGBA<double>         dRGBA;

} /* Close Namespace "gravis" */

#endif
