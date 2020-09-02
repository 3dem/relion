#ifndef __LIBGRAVIS_T_LAB_H__
#define __LIBGRAVIS_T_LAB_H__
/******************************************************************************
**        Title: tLab.h
**  Description: Represents an L*a*b* color tupel.
**
******************************************************************************/

#include <iostream>

namespace gravis
{

  template <class T>
  class tLab
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

      T L, a, b;

      tLab () : L(T(0)), a(T(0)), b(T(0)) { }
      tLab (T L, T a, T b) : L(L), a(a), b(b) { }
      //	tLab (T gray) : (gray), g(gray), b(gray) { }

      void set (T _y, T _cb, T _cr)
      {
        L = _y;
        a = _cb;
        b = _cr;
      }

      //	void add (T _r, T _g, T _b) {
      //		r += _r; g += _g; b += _b;
      //	}

      T intensity () const
      {
        return L;
      }


      /*	bool operator != (const tLab& c) const {
      		return r != c.r || g != c.g || b != c.b;
      	}

      	bool operator == (const tLab& c) const {
      		return r == c.r && g == c.g && b == c.b;
      	}
      */
      tLab& operator += (const tLab& c)
      {
        L += c.L;
        a += c.a;
        b += c.b;
        return *this;
      }
      /*
      	tLab& operator += (const T gray) {
      		r += gray; g += gray; b += gray;
      		return *this;
      	}
      */
      tLab& operator -= (const tLab& c)
      {
        L -= c.L;
        a -= c.a;
        b -= c.b;
        return *this;
      }

      //	tLab& operator -= (const T gray) {
      //		r -= gray; g -= gray; b -= gray;
      //		return *this;
      //	}

      tLab& operator *= (const tLab& c)
      {
        L *= c.L;
        a *= c.a;
        b *= c.b;
        return *this;
      }

      tLab& operator *= (const T factor)
      {
        L *= factor;
        a *= factor;
        b *= factor;
        return *this;
      }
      /*
      	tLab& operator /= (const tLab& c) {
      		r /= c.r; g /= c.g; b /= c.b;
      		return *this;
      	}

      	tLab& operator /= (const T factor) {
      		r /= factor; g /= factor; b /= factor;
      		return *this;
      	}


               * \brief All color components are clamped to [0,1]. This function works inplace.
               *
               * \return self

              tLab& clamp() {
                r = priv::min(priv::max(r, 0), 1);
                g = priv::min(priv::max(g, 0), 1);
                b = priv::min(priv::max(b, 0), 1);
                return *this;
              }

              //! Unary minus
              inline
              tLab operator - () const {
                return tLab<T>(-r, -g, -b);
              };

              //! Addition of a scalar (analog to -=)
              inline
              tLab operator + (const T& c) const { return tLab<T>(r+c, g+c, b+c);  };

              //! Subtraction of a scalar (analog to +=)
              inline
              tLab operator - (const T& c) const { return tLab<T>(r-c, g-c, b-c);  };
             */
      //! Multiplication of a scalar (analog to *=)
      inline
      tLab operator * (const T& c) const
      {
        return tLab<T>(L*c, a*c, b*c);
      };
      /*
              //! Division by a scalar (analog to /=)
              inline
              tLab operator / (const T& c) const { return tLab<T>(r/c, g/c, b/c);  };

      	bool operator == (const tLab& arg) {
      		return ((arg.r == r) && (arg.g == g) && (arg.b == b));
      	}

              const T &operator [](const size_t &i) const { return (&r)[i]; }
              T &operator [](const size_t &i) { return (&r)[i]; }
      */
  };

  template <class T> inline
  tLab<T> operator + (const tLab<T>& c1, const tLab<T>& c2)
  {
    tLab<T> result = c1;
    return (result += c2);
  }

  template <class T> inline
  tLab<T> operator - (const tLab<T>& c1, const tLab<T>& c2)
  {
    tLab<T> result = c1;
    return (result -= c2);
  }

  /*

  template <class T> inline
  tLab<T> operator * (const tLab<T>& c1, const tLab<T>& c2) {
    tLab<T> result(c1.r * c2.r, c1.g * c2.g, c1.b * c2.b);
    return result;
  }
  */
  template <class T> inline
  tLab<T> operator * (const tLab<T>& c, T factor)
  {
    tLab<T> result(c.L * factor, c.a * factor, c.b * factor);
    return result;
  }

  template <class T> inline
  tLab<T> operator * (T factor, const tLab<T>& c)
  {
    tLab<T> result(c.L * factor, c.a * factor, c.b * factor);
    return result;
  }
  /*
  template <class T> inline
  tLab<T> operator / (const tLab<T>& c1, const tLab<T>& c2) {
    tLab<T> result(c1.r / c2.r, c1.g / c2.g, c1.b / c2.b);
    return result;
  }

  template <class T> inline
  tLab<T> operator / (const tLab<T>& c, T factor) {
    tLab<T> result(c.r / factor, c.g / factor, c.b / factor);
    return result;
  }


  template <class T> inline
  bool operator < (const tLab<T>& c1, const tLab<T>& c2) {
    T gray1 = c1.grayValue();
    T gray2 = c2.grayValue();
    return (gray1 < gray2);
  }


  template <class T> inline
  tLab<T> operator ! (const tLab<T>& c) {
    tLab<T> result = tLab<T>::White();
    return (result -= c);
  }

  // Absolute of every color channel
  template <class T> inline
  tLab<T> abs(const tLab<T>& c) {
    return tLab<T>(c.r < T(0) ? -c.r : c.r, c.g < T(0) ? -c.g : c.g, c.b < T(0) ? -c.b : c.b);
  }


  template <class T> inline
  std::ostream& operator << (std::ostream& os, const tLab<T>& c) {
    os << "(" << c.r << " " << c.g << " " << c.b << ")";
    return os;
  }

  template <> inline
  std::ostream& operator << (std::ostream& os, const tLab<unsigned char>& c) {
    os << "(" << (int)c.r << " " << (int)c.g << " " << (int)c.b << ")";
    return os;
  }

  template <class T>
  inline
  T dot (const tLab<T>& v1, const tLab<T>& v2) {
            return (v1.r*v2.r + v1.g*v2.g + v1.b*v2.b);
  }

  */

  //typedef tLab<unsigned char>  bRGB;
  typedef tLab<float>          fLab;
  typedef tLab<double>         dLab;



}

#endif
