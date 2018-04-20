#ifndef __LIBGRAVIS_T_BGR_H__
#define __LIBGRAVIS_T_BGR_H__
/******************************************************************************
**        Title: tBGR.h
**  Description: Represents an BGR color tupel.
**
**       Author:
**
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <iostream>

namespace gravis
{

  template <class T>
  struct tBGR
  {

    T b, g, r;

    typedef T scalar_type;

    tBGR () : b(T(0)), g(T(0)), r(T(0)) { }
    tBGR (T _b, T _g, T _r) : b(_b), g(_g), r(_r) { }
    tBGR (T gray) : b(gray), g(gray), r(gray) { }

    void set (T _b, T _g, T _r)
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

    T grayValue () const
    {
      return (T)(0.30f*r + 0.59f*g + 0.11f*b);
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

    tBGR& operator += (const tBGR& c)
    {
      r += c.r;
      g += c.g;
      b += c.b;
      return *this;
    }

    tBGR& operator += (const T gray)
    {
      r += gray;
      g += gray;
      b += gray;
      return *this;
    }

    tBGR& operator -= (const tBGR& c)
    {
      r -= c.r;
      g -= c.g;
      b -= c.b;
      return *this;
    }

    tBGR& operator -= (const T gray)
    {
      r -= gray;
      g -= gray;
      b -= gray;
      return *this;
    }

    tBGR& operator *= (const tBGR& c)
    {
      r *= c.r;
      g *= c.g;
      b *= c.b;
      return *this;
    }

    tBGR& operator *= (const T factor)
    {
      r *= factor;
      g *= factor;
      b *= factor;
      return *this;
    }

    tBGR& operator /= (const tBGR& c)
    {
      r /= c.r;
      g /= c.g;
      b /= c.b;
      return *this;
    }

    tBGR& operator /= (const T factor)
    {
      r /= factor;
      g /= factor;
      b /= factor;
      return *this;
    }

    //! Unary minus
    inline
    tBGR operator - () const
    {
      return tBGR<T>(-r, -g, -b);
    };

    //! Addition of a scalar (analog to -=)
    inline
    tBGR operator + (const T& c) const
    {
      return tBGR<T>(r+c, g+c, b+c);
    };

    //! Subtraction of a scalar (analog to +=)
    inline
    tBGR operator - (const T& c) const
    {
      return tBGR<T>(r-c, g-c, b-c);
    };

    //! Multiplication of a scalar (analog to *=)
    inline
    tBGR operator * (const T& c) const
    {
      return tBGR<T>(r*c, g*c, b*c);
    };

    //! Division by a scalar (analog to /=)
    inline
    tBGR operator / (const T& c) const
    {
      return tBGR<T>(r/c, g/c, b/c);
    };

    bool operator == (const tBGR& arg)
    {
      return ((arg.r == r) && (arg.g == g) && (arg.b == b));
    }
  };


  template <class T> inline
  tBGR<T> operator + (const tBGR<T>& c1, const tBGR<T>& c2)
  {
    tBGR<T> result = c1;
    return (result += c2);
  }

  template <class T> inline
  tBGR<T> operator - (const tBGR<T>& c1, const tBGR<T>& c2)
  {
    tBGR<T> result = c1;
    return (result -= c2);
  }

  template <class T> inline
  tBGR<T> operator * (const tBGR<T>& c1, const tBGR<T>& c2)
  {
    tBGR<T> result(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b);
    return result;
  }

  template <class T> inline
  tBGR<T> operator * (const tBGR<T>& c, T factor)
  {
    tBGR<T> result(c.r * factor, c.g * factor, c.b * factor);
    return result;
  }

  template <class T> inline
  tBGR<T> operator * (T factor, const tBGR<T>& c)
  {
    tBGR<T> result(c.r * factor, c.g * factor, c.b * factor);
    return result;
  }

  template <class T> inline
  tBGR<T> operator / (const tBGR<T>& c1, const tBGR<T>& c2)
  {
    tBGR<T> result(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b);
    return result;
  }

  template <class T> inline
  tBGR<T> operator / (const tBGR<T>& c, T factor)
  {
    tBGR<T> result(c.r / factor, c.g / factor, c.b / factor);
    return result;
  }


  template <class T> inline
  bool operator < (const tBGR<T>& c1, const tBGR<T>& c2)
  {
    T gray1 = c1.grayValue();
    T gray2 = c2.grayValue();
    return (gray1 < gray2);
  }


  template <class T> inline
  tBGR<T> operator ! (const tBGR<T>& c)
  {
    tBGR<T> result = tBGR<T>::White();
    return (result -= c);
  }

  // Absolute of every color channel
  template <class T> inline
  tBGR<T> abs(const tBGR<T>& c)
  {
    return tBGR<T>(c.r < T(0) ? -c.r : c.r, c.g < T(0) ? -c.g : c.g, c.b < T(0) ? -c.b : c.b);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tBGR<T>& c)
  {
    os << "(" << c.r << " " << c.g << " " << c.b << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tBGR<unsigned char>& c)
  {
    os << "(" << (int)c.r << " " << (int)c.g << " " << (int)c.b << ")";
    return os;
  }

  typedef tBGR<char>           cBGR;
  typedef tBGR<unsigned char>  bBGR;
  typedef tBGR<float>          fBGR;
  typedef tBGR<double>         dBGR;

} /* Close Namespace "gravis" */

#endif
