#ifndef __LIBGRAVIS_T_YCB_CR_H__
#define __LIBGRAVIS_T_YCB_CR_H__
/******************************************************************************
**        Title: tYCbCr.h
**  Description: Represents an CIE Y/Cb/Cr color tupel.
**
******************************************************************************/

#include <iostream>

namespace gravis
{

  template <class T>
  class tYCbCr
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

      T y, cb, cr;

      tYCbCr () : y(T(0)), cb(T(0)), cr(T(0)) { }
      tYCbCr (T y, T cb, T cr) : y(y), cb(cb), cr(cr) { }
      //	tYCbCr (T gray) : (gray), g(gray), b(gray) { }

      void set (T _y, T _cb, T _cr)
      {
        y = _y;
        cb = _cb;
        cr = _cr;
      }

      //	void add (T _r, T _g, T _b) {
      //		r += _r; g += _g; b += _b;
      //	}

      T intensity () const
      {
        return y();
      }


      /*	bool operator != (const tYCbCr& c) const {
      		return r != c.r || g != c.g || b != c.b;
      	}

      	bool operator == (const tYCbCr& c) const {
      		return r == c.r && g == c.g && b == c.b;
      	}
      */
      tYCbCr& operator += (const tYCbCr& c)
      {
        y += c.y;
        cb += c.cb;
        cr += c.cr;
        return *this;
      }
      /*
      	tYCbCr& operator += (const T gray) {
      		r += gray; g += gray; b += gray;
      		return *this;
      	}
      */
      tYCbCr& operator -= (const tYCbCr& c)
      {
        y -= c.y;
        cb -= c.cb;
        cr -= c.cr;
        return *this;
      }

      //	tYCbCr& operator -= (const T gray) {
      //		r -= gray; g -= gray; b -= gray;
      //		return *this;
      //	}

      tYCbCr& operator *= (const tYCbCr& c)
      {
        y *= c.y;
        cb *= c.cb;
        cr *= c.cr;
        return *this;
      }

      tYCbCr& operator *= (const T factor)
      {
        y *= factor;
        cb *= factor;
        cr *= factor;
        return *this;
      }
      /*
      	tYCbCr& operator /= (const tYCbCr& c) {
      		r /= c.r; g /= c.g; b /= c.b;
      		return *this;
      	}

      	tYCbCr& operator /= (const T factor) {
      		r /= factor; g /= factor; b /= factor;
      		return *this;
      	}


               * \brief All color components are clamped to [0,1]. This function works inplace.
               *
               * \return self

              tYCbCr& clamp() {
                r = std::min(std::max(r, 0), 1);
                g = std::min(std::max(g, 0), 1);
                b = std::min(std::max(b, 0), 1);
                return *this;
              }

              //! Unary minus
              inline
              tYCbCr operator - () const {
                return tYCbCr<T>(-r, -g, -b);
              };

              //! Addition of a scalar (analog to -=)
              inline
              tYCbCr operator + (const T& c) const { return tYCbCr<T>(r+c, g+c, b+c);  };

              //! Subtraction of a scalar (analog to +=)
              inline
              tYCbCr operator - (const T& c) const { return tYCbCr<T>(r-c, g-c, b-c);  };
             */
      //! Multiplication of a scalar (analog to *=)
      inline
      tYCbCr operator * (const T& c) const
      {
        return tYCbCr<T>(y*c, cb*c, cr*c);
      };
      /*
              //! Division by a scalar (analog to /=)
              inline
              tYCbCr operator / (const T& c) const { return tYCbCr<T>(r/c, g/c, b/c);  };

      	bool operator == (const tYCbCr& arg) {
      		return ((arg.r == r) && (arg.g == g) && (arg.b == b));
      	}

              const T &operator [](const size_t &i) const { return (&r)[i]; }
              T &operator [](const size_t &i) { return (&r)[i]; }
      */
  };

  template <class T> inline
  tYCbCr<T> operator + (const tYCbCr<T>& c1, const tYCbCr<T>& c2)
  {
    tYCbCr<T> result = c1;
    return (result += c2);
  }

  template <class T> inline
  tYCbCr<T> operator - (const tYCbCr<T>& c1, const tYCbCr<T>& c2)
  {
    tYCbCr<T> result = c1;
    return (result -= c2);
  }

  /*

  template <class T> inline
  tYCbCr<T> operator * (const tYCbCr<T>& c1, const tYCbCr<T>& c2) {
    tYCbCr<T> result(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b);
    return result;
  }
  */
  template <class T> inline
  tYCbCr<T> operator * (const tYCbCr<T>& c, T factor)
  {
    tYCbCr<T> result(c.y * factor, c.cb * factor, c.cr * factor);
    return result;
  }

  template <class T> inline
  tYCbCr<T> operator * (T factor, const tYCbCr<T>& c)
  {
    tYCbCr<T> result(c.y * factor, c.cb * factor, c.cr * factor);
    return result;
  }
  /*
  template <class T> inline
  tYCbCr<T> operator / (const tYCbCr<T>& c1, const tYCbCr<T>& c2) {
    tYCbCr<T> result(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b);
    return result;
  }

  template <class T> inline
  tYCbCr<T> operator / (const tYCbCr<T>& c, T factor) {
    tYCbCr<T> result(c.r / factor, c.g / factor, c.b / factor);
    return result;
  }


  template <class T> inline
  bool operator < (const tYCbCr<T>& c1, const tYCbCr<T>& c2) {
    T gray1 = c1.grayValue();
    T gray2 = c2.grayValue();
    return (gray1 < gray2);
  }


  template <class T> inline
  tYCbCr<T> operator ! (const tYCbCr<T>& c) {
    tYCbCr<T> result = tYCbCr<T>::White();
    return (result -= c);
  }

  // Absolute of every color channel
  template <class T> inline
  tYCbCr<T> abs(const tYCbCr<T>& c) {
    return tYCbCr<T>(c.r < T(0) ? -c.r : c.r, c.g < T(0) ? -c.g : c.g, c.b < T(0) ? -c.b : c.b);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tYCbCr<T>& c) {
    os << "(" << c.r << " " << c.g << " " << c.b << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tYCbCr<unsigned char>& c) {
    os << "(" << (int)c.r << " " << (int)c.g << " " << (int)c.b << ")";
    return os;
  }

  template <class T>
  inline
  T dot (const tYCbCr<T>& v1, const tYCbCr<T>& v2) {
            return (v1.r*v2.r + v1.g*v2.g + v1.b*v2.b);
  }

  */

  //typedef tYCbCr<unsigned char>  bRGB;
  typedef tYCbCr<float>          fYCbCr;
  typedef tYCbCr<double>         dYCbCr;



}

#endif
