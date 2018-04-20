#ifndef __LIBGRAVIS_T_RGB_H__
#define __LIBGRAVIS_T_RGB_H__
/******************************************************************************
**        Title: tRGB.h
**  Description: Represents an RGB color tupel.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Brian Amberg, 2005-2006
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <iostream>
#include <stdexcept>
#include <sstream>

namespace gravis
{

  template <class T> class tRGBA;

  template <class T>
  class tRGB
  {
      /*!
       * Private helper functions, wrapped into an additional struct in case that we want to use the names
       **/
      struct priv
      {
        static inline const T& min(const T& a, const T& b)
        {
          return a<b ? a : b;
        }
        static inline const T& max(const T& a, const T& b)
        {
          return a>b ? a : b;
        }
      };
    public:

      typedef T scalar_type;

      T r, g, b;

      tRGB () : r(T(0)), g(T(0)), b(T(0)) { }
      tRGB (T _r, T _g, T _b) : r(_r), g(_g), b(_b) { }
      tRGB (T gray) : r(gray), g(gray), b(gray) { }
      explicit tRGB (const tRGBA<T>& c) : r(c.r), g(c.g), b(c.b) {}

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

      void add (T _r, T _g, T _b)
      {
        r += _r;
        g += _g;
        b += _b;
      }

      void add (T gray)
      {
        r += gray;
        g += gray;
        b += gray;
      }

      /**
       * Deprecated, use intensity() instead
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
       * Return minimum of the tupel.
       *
       * TODO: Is this really necessary in here. It could be a utility function.
       **/
      T const& minValue () const
      {
        return std::min(std::min(r, g), b);
      }

      /*!
       * Return maximum of the tupel.
       *
       * TODO: Is this really necessary in here. It could be a utility function.
       **/
      T const& maxValue () const
      {
        return std::max(std::max(r, g), b);
      }

      bool operator != (const tRGB& c) const
      {
        return r != c.r || g != c.g || b != c.b;
      }

      bool operator == (const tRGB& c) const
      {
        return r == c.r && g == c.g && b == c.b;
      }

      tRGB& operator += (const tRGB& c)
      {
        r += c.r;
        g += c.g;
        b += c.b;
        return *this;
      }

      tRGB& operator += (const T gray)
      {
        r += gray;
        g += gray;
        b += gray;
        return *this;
      }

      tRGB& operator -= (const tRGB& c)
      {
        r -= c.r;
        g -= c.g;
        b -= c.b;
        return *this;
      }

      tRGB& operator -= (const T gray)
      {
        r -= gray;
        g -= gray;
        b -= gray;
        return *this;
      }

      tRGB& operator *= (const tRGB& c)
      {
        r *= c.r;
        g *= c.g;
        b *= c.b;
        return *this;
      }

      tRGB& operator *= (const T factor)
      {
        r *= factor;
        g *= factor;
        b *= factor;
        return *this;
      }

      tRGB& operator /= (const tRGB& c)
      {
        r /= c.r;
        g /= c.g;
        b /= c.b;
        return *this;
      }

      tRGB& operator /= (const T factor)
      {
        r /= factor;
        g /= factor;
        b /= factor;
        return *this;
      }

      /*!
       * \brief All color components are clamped to [0,1]. This function works inplace.
       *
       * \return self
       */
      tRGB& clamp()
      {
        r = std::min(std::max(r, T(0)), T(1));
        g = std::min(std::max(g, T(0)), T(1));
        b = std::min(std::max(b, T(0)), T(1));
        return *this;
      }

      //! Unary minus
      inline
      tRGB operator - () const
      {
        return tRGB<T>(-r, -g, -b);
      };

      //! Addition of a scalar (analog to -=)
      inline
      tRGB operator + (const T& c) const
      {
        return tRGB<T>(r+c, g+c, b+c);
      };

      //! Subtraction of a scalar (analog to +=)
      inline
      tRGB operator - (const T& c) const
      {
        return tRGB<T>(r-c, g-c, b-c);
      };

      //! Multiplication of a scalar (analog to *=)
      inline
      tRGB operator * (const T& c) const
      {
        return tRGB<T>(r*c, g*c, b*c);
      };

      //! Division by a scalar (analog to /=)
      inline
      tRGB operator / (const T& c) const
      {
        return tRGB<T>(r/c, g/c, b/c);
      };

      bool operator == (const tRGB& arg)
      {
        return ((arg.r == r) && (arg.g == g) && (arg.b == b));
      }

      const T& operator [](const size_t& i) const
      {
        return (&r)[i];
      }
      T& operator [](const size_t& i)
      {
        return (&r)[i];
      }
  };


  template <class T> inline
  tRGB<T> operator + (const tRGB<T>& c1, const tRGB<T>& c2)
  {
    tRGB<T> result = c1;
    return (result += c2);
  }

  template <class T> inline
  tRGB<T> operator - (const tRGB<T>& c1, const tRGB<T>& c2)
  {
    tRGB<T> result = c1;
    return (result -= c2);
  }

  template <class T> inline
  tRGB<T> operator * (const tRGB<T>& c1, const tRGB<T>& c2)
  {
    tRGB<T> result(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b);
    return result;
  }

  template <class T> inline
  tRGB<T> operator * (const tRGB<T>& c, T factor)
  {
    tRGB<T> result(c.r * factor, c.g * factor, c.b * factor);
    return result;
  }

  template <class T> inline
  tRGB<T> operator * (T factor, const tRGB<T>& c)
  {
    tRGB<T> result(c.r * factor, c.g * factor, c.b * factor);
    return result;
  }

  template <class T> inline
  tRGB<T> operator / (const tRGB<T>& c1, const tRGB<T>& c2)
  {
    tRGB<T> result(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b);
    return result;
  }

  template <class T> inline
  tRGB<T> operator / (const tRGB<T>& c, T factor)
  {
    tRGB<T> result(c.r / factor, c.g / factor, c.b / factor);
    return result;
  }


  template <class T> inline
  bool operator < (const tRGB<T>& c1, const tRGB<T>& c2)
  {
    T gray1 = c1.grayValue();
    T gray2 = c2.grayValue();
    return (gray1 < gray2);
  }


  template <class T> inline
  tRGB<T> operator ! (const tRGB<T>& c)
  {
    tRGB<T> result = tRGB<T>::White();
    return (result -= c);
  }

  // Absolute of every color channel
  template <class T> inline
  tRGB<T> abs(const tRGB<T>& c)
  {
    return tRGB<T>(c.r < T(0) ? -c.r : c.r, c.g < T(0) ? -c.g : c.g, c.b < T(0) ? -c.b : c.b);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tRGB<T>& c)
  {
    os << "(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tRGB<unsigned char>& c)
  {
    os << "(" << (int)c.r << ", " << (int)c.g << ", " << (int)c.b << ")";
    return os;
  }

    // Inverse of operator<<
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, tRGB<T>& arg)
  {
    char c = ' ';
    is >> c;
    if (is.eof())
      return is;
    if (c != '(')
      throw std::runtime_error("tRGB should start with an opening (");
    std::string values;
    int v = 0;
    while ((is.get(c)) && (c != ')'))
    {
      if (c == ',')
      {
        v++;
        if (v >= 3)
          throw std::runtime_error("tRGB contains more than three elements");
        values.push_back(' ');
      }
      else
        values.push_back(c);
    }
    if (c != ')')
    {
      throw std::runtime_error("tRGB should end with a )");
    }
    if ( v < 2 )
    {
      throw std::runtime_error("tRGB has not enough color values");
    }
      
    std::stringstream valueReader(values);
    valueReader >> arg.r >> arg.g >> arg.b;
    return is;
  }

  
  
  template <class T>
  inline
  T dot (const tRGB<T>& v1, const tRGB<T>& v2)
  {
    return (v1.r*v2.r + v1.g*v2.g + v1.b*v2.b);
  }


  typedef tRGB<unsigned char>  bRGB;
  typedef tRGB<float>          fRGB;
  typedef tRGB<double>         dRGB;

} /* Close Namespace "gravis" */

#endif
