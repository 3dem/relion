#ifndef __LIBGRAVIS_T_MATRIX_H__
#define __LIBGRAVIS_T_MATRIX_H__
/******************************************************************************
**	  Title: matrix.h
**  Description: Templated fixed size dense matrices, which are a
**               complement to the fixed size t{2,3,4}{Vector,Matrix} classes.
**
**	 Author: Brian Amberg, 2007
**		 Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "StringFormat.h"
#include "Exception.h"

#ifdef DEBUG
#define checkAccess1( i,     h     ) { if(!( (i)<(h)                      )) GRAVIS_THROW2(gravis::Exception, "Access out of bounds " #i "<" #h);                                 }
#define checkAccess2( i,j,   h,w   ) { if(!( (i)<(h) && (j)<(w)           )) GRAVIS_THROW2(gravis::Exception, "Access out of bounds " #i "<" #h " && "#j"<" #w);                  }
#define checkAccess3( i,j,k, h,w,d ) { if(!( (i)<(h) && (j)<(w) && (k)<(d))) GRAVIS_THROW2(gravis::Exception, "Access out of bounds " #i "<" #h " && "#j"<" #w " && " #k "<" #d); }
#else
#define checkAccess1( i,     h     ) { }
#define checkAccess2( i,j,   h,w   ) { }
#define checkAccess3( i,j,k, h,w,d ) { }
#endif


namespace gravis
{

  template <class T, size_t h, size_t w> class tMatrix;

  namespace tMatrixPrivateConstructorTrick
  {
    struct CheckIfRightSize
    {
      template <class T> static void has_2_elements( const tMatrix<T, 2, 1> &m ) {}
      template <class T> static void has_2_elements( const tMatrix<T, 1, 2> &m ) {}

      template <class T> static void has_3_elements( const tMatrix<T, 3, 1> &m ) {}
      template <class T> static void has_3_elements( const tMatrix<T, 1, 3> &m ) {}

      template <class T> static void has_4_elements( const tMatrix<T, 4, 1> &m ) {}
      template <class T> static void has_4_elements( const tMatrix<T, 1, 4> &m ) {}
      template <class T> static void has_4_elements( const tMatrix<T, 2, 2> &m ) {}

      template <class T> static void has_9_elements( const tMatrix<T, 3, 3> &m ) {}

      template <class T> static void has_16_elements( const tMatrix<T, 4, 4> &m ) {}
    };
  }


  /**
   * Small Matrix Of Arbitrary size held completely in memory in consecutive positions.
   * The data is in row major order.
   **/
  template <class T, size_t h, size_t w>
  class tMatrix
  {
    public:
      typedef T scalar;
      T data[h* w];

      /**
       * The data is not initialized
       **/
      tMatrix() {};
      /**
       * Copy constructor
       **/
      tMatrix(const tMatrix<T, h, w> &o)
      {
        memcpy( data, o.data, h*w*sizeof(T) );
      }
      /**
       * Fill with copies of v values
       **/
      explicit tMatrix(const T& v)
      {
        fill(v);
      }

      /**
       * Copy data from another matrix
       **/
      tMatrix& operator=(const tMatrix& o)
      {
        memcpy( data, o.data, h*w*sizeof(T) );
        return *this;
      }

      /**
       * Inplace negation
       **/
      inline
      void negate()
      {
        tMatrix<T, h, w> &m = *this;
        for (size_t i=0; i<m.size(); ++i) m[i] = -m[i];
      }

      /**
       * Special constructors for certain sizes of small matrices. It's
       * constructed to fail compiling when a wrong size is used.
       *
       * Constructor for a 2x1 or 1x2 vector
       **/
      tMatrix(const T& e0, const T& e1)
      {
        tMatrixPrivateConstructorTrick::CheckIfRightSize::has_2_elements(*this);
        data[0] = e0;
        data[1] = e1;
      }

      /**
       * Special constructors for certain sizes of small matrices. It's
       * constructed to fail compiling when a wrong size is used.
       *
       * Constructor for a 3x1 or 1x3 vector
       **/
      tMatrix(const T& e0, const T& e1, const T& e2)
      {
        tMatrixPrivateConstructorTrick::CheckIfRightSize::has_3_elements(*this);
        data[0] = e0;
        data[1] = e1;
        data[2] = e2;
      }

      /**
       * Special constructors for certain sizes of small matrices. It's
       * constructed to fail compiling when a wrong size is used.
       *
       * Constructor for a 4x1 or 1x4 vector, or a 2x2 matrix. (Beware, the elements are row major, so they seem to be transposed.
       **/
      tMatrix(const T& e0, const T& e1, const T& e2, const T& e3)
      {
        tMatrixPrivateConstructorTrick::CheckIfRightSize::has_4_elements(*this);
        data[0] = e0;
        data[1] = e1;
        data[2] = e2;
        data[3] = e3;
      }

      /**
       * Special constructors for certain sizes of small matrices. It's
       * constructed to fail compiling when a wrong size is used.
       *
       * Constructor for a 3x3 matrix. (Beware, the elements are row major, so they seem to be transposed.
       **/
      tMatrix(
        const T& e0, const T& e1, const T& e2,
        const T& e3, const T& e4, const T& e5,
        const T& e6, const T& e7, const T& e8
      )
      {
        tMatrixPrivateConstructorTrick::CheckIfRightSize::has_9_elements(*this);
        data[0] = e0;
        data[1] = e1;
        data[2] = e2;
        data[3] = e3;
        data[4] = e4;
        data[5] = e5;
        data[6] = e6;
        data[7] = e7;
        data[7] = e8;
      }

      /**
       * Special constructors for certain sizes of small matrices. It's
       * constructed to fail compiling when a wrong size is used.
       *
       * Constructor for a 4x4 matrix. (Beware, the elements are row major, so they seem to be transposed.
       **/
      tMatrix(
        const T& e0,  const T& e1,  const T& e2,  const T& e3,
        const T& e4,  const T& e5,  const T& e6,  const T& e7,
        const T& e8,  const T& e9,  const T& e10, const T& e11,
        const T& e12, const T& e13, const T& e14, const T& e15
      )
      {
        tMatrixPrivateConstructorTrick::CheckIfRightSize::has_16_elements(*this);
        data[0]  = e0;
        data[1]  = e1;
        data[2]  = e2;
        data[3] = e3;
        data[4]  = e4;
        data[5]  = e5;
        data[6]  = e6;
        data[7] = e7;
        data[8]  = e8;
        data[9]  = e9;
        data[10] = e10;
        data[11] = e11;
        data[12] = e12;
        data[13] = e13;
        data[14] = e14;
        data[15] = e15;
      }

      /**
       * Number of elements
       **/
      inline size_t size() const
      {
        return h*w;
      };

      /**
       * Access is checked when compiled with DEBUG flag
       * The data is in row major order
       **/
      inline const T& operator[](size_t i) const
      {
        checkAccess1(i, h*w);
        return data[i];
      };
      /**
       * Access is checked when compiled with DEBUG flag
       * The data is in row major order
       **/
      inline       T& operator[](size_t i)
      {
        checkAccess1(i, h*w);
        return data[i];
      };

      /**
       * Access is checked when compiled with DEBUG flag
       * The data is in row major order
       **/
      inline const T& operator()(size_t r, size_t c) const
      {
        checkAccess2(r, c, h, w);
        return data[r+c*h];
      };
      /**
       * Access is checked when compiled with DEBUG flag
       * The data is in row major order
       **/
      inline       T& operator()(size_t r, size_t c)
      {
        checkAccess2(r, c, h, w);
        return data[r+c*h];
      };

      /**
       * Negate a matrix
       **/
      inline tMatrix operator - (void)
      {
        tMatrix t;
        for (size_t i=0; i<size(); ++i) t[i] = -(*this)[i];
        return t;
      }

      /**
       * Dot product
       **/
      inline T dot(const tMatrix& o) const
      {
        T r((*this)[0]*o[0]);
        for (size_t i=1; i<size(); ++i)
          r += (*this)[i]*o[i];
        return r;
      };

      /**
       * Fill with copies
       **/
      inline void fill(const T& v)
      {
        for (size_t i=0; i<size(); ++i) (*this)[i] = v;
      }
      /**
       * Fill with ones
       **/
      inline void ones()
      {
        fill(T(1));
      }
      /**
       * Fill with zeros
       **/
      inline void zeros()
      {
        fill(T(0));
      }
      /**
       * Set to identity
       **/
      inline void identity()
      {
        fill(0);
        if (w<h)
          for (size_t i=0; i<w; ++i) (*this)(i,i) = 1;
        else
          for (size_t i=0; i<h; ++i) (*this)(i,i) = 1;
      }

      /**
       * Matrix Scalar Multiplication
       **/
      tMatrix  operator* (const T& o) const
      {
        tMatrix r;
        for (size_t i=0; i<size(); ++i) r[i] = (*this)[i] * o;
        return r;
      }
      /**
       * Matrix Scalar Division
       **/
      tMatrix  operator/ (const T& o) const
      {
        tMatrix r;
        for (size_t i=0; i<size(); ++i) r[i] = (*this)[i] / o;
        return r;
      }
      /**
       * Matrix Scalar Multiplication
       **/
      tMatrix& operator*=(const T& o)
      {
        for (size_t i=0; i<size(); ++i) (*this)[i] *= o;
        return *this;
      }
      /**
       * Matrix Scalar Division
       **/
      tMatrix& operator/=(const T& o)
      {
        for (size_t i=0; i<size(); ++i) (*this)[i] /= o;
        return *this;
      }
      /**
       * Matrix-Matrix Multiplication
       **/
      template <size_t w2>
      tMatrix<T,h,w2> operator*(const tMatrix<T, w,w2> &right) const
      {
        tMatrix<T,h,w2> out(0);
        const tMatrix& self(*this);
        for (size_t j=0; j<w2; ++j)
          for (size_t k=0; k<w; ++k)
            for (size_t i=0; i<h; ++i)
              out(i,j) += self(i, k) * right(k, j);
        return out;
      }
      /**
       * Matrix-Scalar Addition
       **/
      tMatrix& operator+=(const T& o)
      {
        for (size_t i=0; i<size(); ++i) (*this)[i] += o;
        return *this;
      }
      /**
       * Matrix-Scalar Addition
       **/
      tMatrix operator+(const T& o) const
      {
        tMatrix r;
        for (size_t i=0; i<size(); ++i) r[i] = (*this)[i] + o;
        return r;
      }
      /**
       * Matrix-Matrix Addition
       **/
      tMatrix& operator+=(const tMatrix& o)
      {
        for (size_t i=0; i<size(); ++i) (*this)[i] += o[i];
        return *this;
      }
      /**
       * Matrix-Matrix Addition
       **/
      tMatrix operator+(const tMatrix& o) const
      {
        tMatrix r;
        for (size_t i=0; i<size(); ++i) r[i] = (*this)[i] + o[i];
        return r;
      }
      /**
       * Matrix-Scalar Subtraction
       **/
      tMatrix& operator-=(const T& o)
      {
        for (size_t i=0; i<size(); ++i) (*this)[i] -= o;
        return *this;
      }
      /**
       * Matrix-Scalar Subtraction
       **/
      tMatrix operator-(const T& o) const
      {
        tMatrix r;
        for (size_t i=0; i<size(); ++i) r[i] = (*this)[i] - o;
        return r;
      }
      /**
       * Matrix-Matrix Subtraction
       **/
      tMatrix& operator-=(const tMatrix& o)
      {
        sub(*this, o);
        return *this;
      }
      /**
       * Matrix-Matrix Subtraction
       **/
      tMatrix operator-(const tMatrix& o) const
      {
        tMatrix r;
        for (size_t i=0; i<size(); ++i) r[i] = (*this)[i] - o[i];
        return r;
      }

      /**
       * Interpret matrix as scalar of underlying type, returns the first entry
       **/
      // operator       T()       { return (*this)[0]; };
      /**
       * L2 Norm
       **/
      T normL2sqr() const
      {
        if (size()==0) return 0;
        T result = (*this)[0]*(*this)[0];
        for (size_t i=1; i<size(); ++i)
          result += (*this)[i]*(*this)[i];
        return result;
      }
      /**
       * L2 Norm
       **/
      T normL2() const
      {
        return sqrt(normL2sqr());
      }

      /**
       * Convenience function to clear a matrix
       **/
      inline void clear()
      {
        if (size()>0) memset( &data[0], 0, sizeof(data[0])*size());
      }
      /**
       * Convenience function to clamp all elements of
       **/
      inline void clamp(const T& min, const T& max)
      {
        for (size_t i=0; i<size(); ++i)
        {
          if (data[i]<min) data[i] = min;
          if (max<data[i]) data[i] = max;
        }
      }

  };

  namespace matrix
  {
    /**
     * 3-Vector Cross Product
     **/
    template <class T>
    static inline
    void cross(tMatrix<T,3,1> &result, const tMatrix<T,3,1>& a, const tMatrix<T,3,1>& b)
    {
      result[0] = a[1]*b[2] - a[2]*b[1];
      result[1] = a[2]*b[0] - a[0]*b[2];
      result[2] = a[0]*b[1] - a[1]*b[0];
    }
  }
  /**
   * Matrix Scalar Addition
   **/
  template <class T, size_t h, size_t w>
  inline
  tMatrix<T,h,w> operator-(const T& o, const tMatrix<T, h, w> &self)
  {
    tMatrix<T,h,w> r;
    for (size_t i=0; i<self.size(); ++i) r[i] = o-self[i];
    return r;
  }
  /**
   * Matrix Scalar Addition
   **/
  template <class T, size_t h, size_t w>
  inline
  tMatrix<T,h,w> operator+(const T& o, const tMatrix<T, h, w> &self)
  {
    tMatrix<T,h,w> r;
    for (size_t i=0; i<self.size(); ++i) r[i] = o+self[i];
    return r;
  }
  /**
   * Matrix Scalar Multiplication
   **/
  template <class T, size_t h, size_t w>
  inline
  tMatrix<T,h,w> operator*(const T& o, const tMatrix<T, h, w> &self)
  {
    tMatrix<T,h,w> r;
    for (size_t i=0; i<self.size(); ++i) r[i] = o*self[i];
    return r;
  }

  // Operations on matrices

  /**
   * out += left * right
   **/
  template <class T, size_t h, size_t w, size_t w2>
  inline static
  void addmult(tMatrix<T, h, w2> &out, const tMatrix<T, h, w> &left, const tMatrix<T, w, w2> &right)
  {
    for (size_t j=0; j<w2; ++j)
      for (size_t k=0; k<w; ++k)
        for (size_t i=0; i<h; ++i)
          out(i,j) += left(i, k) * right(k, j);
  }

  /**
   * out -= left * right
   **/
  template <class T, size_t h, size_t w, size_t w2>
  inline static
  void submult(tMatrix<T, h, w2> &out, const tMatrix<T, h, w> &left, const tMatrix<T, w, w2> &right)
  {
    for (size_t j=0; j<w2; ++j)
      for (size_t k=0; k<w; ++k)
        for (size_t i=0; i<h; ++i)
          out(i,j) -= left(i, k) * right(k, j);
  }

  /**
   * Inplace negation
   **/
  template <class T, size_t h, size_t w>
  inline static
  void negate(tMatrix<T, h, w> &m)
  {
    for (size_t i=0; i<m.size(); ++i) m[i] = -m[i];
  }

  /**
   * Inplace scalar matrix multiplication left
   **/
  template <class T, size_t h, size_t w>
  inline static
  void mult(const T& scalar, tMatrix<T, h, w> &m)
  {
    for (size_t i=0; i<m.size(); ++i)
      m[i] = scalar * m[i];
  }

  /**
   * Inplace matrix scalar multiplication
   **/
  template <class T, size_t h, size_t w>
  inline static
  void mult(tMatrix<T, h, w> &m, const T& scalar)
  {
    for (size_t i=0; i<m.size(); ++i)
      m[i] *= scalar;
  }

  /**
   * Matrix scalar multiplication
   **/
  template <class T, size_t h, size_t w>
  inline static
  void mult(tMatrix<T, h, w> &out, const tMatrix<T, h, w> &m, const T& scalar)
  {
    for (size_t i=0; i<m.size(); ++i)
      out[i] = m[i] * scalar;
  }

  /**
   * Scalar Matrix multiplication
   **/
  template <class T, size_t h, size_t w>
  inline static
  void mult(tMatrix<T, h, w> &out, const T& scalar, const tMatrix<T, h, w> &m)
  {
    for (size_t i=0; i<m.size(); ++i)
      out[i] = scalar * m[i];
  }

  /**
   * Inplace Matrix addition
   **/
  template <class T, size_t h, size_t w>
  inline
  void add(tMatrix<T, h, w> &self, const tMatrix<T, h, w> &right)
  {
    for (size_t i=0; i<h*w; ++i)
      self[i]+=right[i];
  }

  /**
   * Matrix addition
   **/
  template <class T, size_t h, size_t w>
  inline
  void add(tMatrix<T, h, w> &out, const tMatrix<T, h, w> &self, const tMatrix<T, h, w> &right)
  {
    for (size_t i=0; i<h*w; ++i)
      out[i] = self[i]+right[i];
  }

  /**
   * Inplace matrix subtraction
   **/
  template <class T, size_t h, size_t w>
  inline
  void sub(tMatrix<T, h, w> &self, const tMatrix<T, h, w> &right)
  {
    for (size_t i=0; i<h*w; ++i)
      self[i]-=right[i];
  }

  /**
   * Matrix subtraction
   **/
  template <class T, size_t h, size_t w>
  inline
  void sub(tMatrix<T, h, w> &out, const tMatrix<T, h, w> &self, const tMatrix<T, h, w> &right)
  {
    for (size_t i=0; i<h*w; ++i)
      out[i] = self[i]-right[i];
  }

  /**
   * Matrix multipliciation
   **/
  template <class T, size_t h, size_t w, size_t w2>
  inline static
  void mult(tMatrix<T,h,w2> &out, const tMatrix<T,h,w> &self, const tMatrix<T, w,w2> &right)
  {
    out.zeros();
    for (size_t i=0; i<h; ++i)
      for (size_t j=0; j<w2; ++j)
        for (size_t k=0; k<w; ++k)
          out(i,j) += self(i, k) * right(k, j);
  }

  /**
   * Vector multiplication (i.e. dot product)
   **/
  template <class T, size_t w>
  inline static
  void mult(T& out, const tMatrix<T,1,w> &self, const tMatrix<T, w,1> &right)
  {
    out = self[0] * right[0];
    for (size_t i=1; i<w; ++i)
      out += self[i] * right[i];
  }

  /**
   * Inplace matrix multiplication for quadratic matrices
   **/
  template <class T, size_t h>
  inline static
  void mult(tMatrix<T,h,h> &self, const tMatrix<T, h,h> &right)
  {
    tMatrix<T,h,h> tmp;
    mult(tmp, self, right);
    self = tmp;
  }


  /** Convenience Constructors **/
  template <class T> inline static tMatrix<T, 1, 1> tVector1(const T& a)
  {
    tMatrix<T, 1, 1> v;
    v[0]=a;
    return v;
  }
  /** Convenience Constructors **/
  template <class T> inline static tMatrix<T, 2, 1> tVector2(const T& a, const T& b)
  {
    tMatrix<T, 2, 1> v;
    v[0]=a;
    v[1]=b;
    return v;
  }
  /** Convenience Constructors **/
  template <class T> inline static tMatrix<T, 3, 1> tVector3(const T& a, const T& b, const T& c)
  {
    tMatrix<T, 3, 1> v;
    v[0]=a;
    v[1]=b;
    v[2]=c;
    return v;
  }
  /** Convenience Constructors **/
  template <class T> inline static tMatrix<T, 4, 1> tVector4(const T& a, const T& b, const T& c, const T& d)
  {
    tMatrix<T, 4, 1> v;
    v[0]=a;
    v[1]=b;
    v[2]=c;
    v[3]=d;
    return v;
  }

  /** Convenience Constructors **/
  template <class T> inline static tMatrix<T, 3, 3> tMatrix3(
    const T& a, const T& b, const T& c,
    const T& d, const T& e, const T& f,
    const T& g, const T& i, const T& h)
  {
    tMatrix<T, 3, 3> m;
    m(0,0)=a;
    m(0,1)=b;
    m(0,2)=c;
    m(1,0)=d;
    m(1,1)=e;
    m(1,2)=f;
    m(2,0)=g;
    m(2,1)=i;
    m(2,2)=h;
    return m;
  }


  /**
   * Write fixed size matrices to a stream
   **/
  template <class T, size_t h, size_t w>
  inline
  std::ostream& operator<< (std::ostream& os, const tMatrix<T, h, w>& arg)
  {
    if ((h>1) && (w>1))
    {
      os << "Matrix: " << h << "x" << w << std::endl;
      for (size_t i=0; i<h; ++i)
      {
        if (i==0)
        {
          os << "/ ";
        }
        else if (i==h-1)
        {
          os << " \\";
        }
        else
        {
          os << " |";
        }
        for (size_t j=0; j<w; ++j) os << std::setw(8) << arg(i,j) << " ";
        if (i==0)
        {
          os << "\\ ";
        }
        else if (i==h-1)
        {
          os << " /";
        }
        else
        {
          os << " |";
        }
        os << "\n";
      }
    }
    else if (w==1 && h>1)
    {
      os << "[ ";
      for (size_t j=0; j<h; ++j)
        os << std::setw(8) << arg[j] << " ";
      os << " ]^T";
    }
    else
    {
      os << "[ ";
      for (size_t j=0; j<w; ++j)
        os << std::setw(8) << arg[j] << " ";
      os << " ]";
    }
    return os;
  }


  /**
   * Read fixed size matrices from a stream
   **/
  template <class T, size_t h, size_t w>
  inline
  std::istream& operator>> (std::istream& is, tMatrix<T, h, w>& arg)
  {
    std::string t;
    if ((h>1) && (w>1))
    {
      is >> t >> t;
      for (size_t i=0; i<h; ++i)
      {
        is >> t;
        for (size_t j=0; j<w; ++j) is >> arg(i,j);
        is >> t;
      }
    }
    else if (w==1 && h>1)
    {
      is >> t;
      if (t != "[") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A vector should start with [", t);
      for (size_t j=0; j<h; ++j)
        is >> arg[j];
      is >> t;
      if (t != "]^T") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A column vector should end with ]^T", t);
    }
    else
    {
      is >> t;
      if (t != "[") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A vector should start with [", t);
      for (size_t j=0; j<h; ++j)
        is >> arg[j];
      is >> t;
      if (t != "]") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A row vector should end with ]", t);
    }
    return is;
  }

}

#endif
