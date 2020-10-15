#ifndef __LIBGRAVIS_T_VAR_MATRIX_H__
#define __LIBGRAVIS_T_VAR_MATRIX_H__
/******************************************************************************
**	  Title: matrix.h
**  Description: Templated variable size dense matrices, with a blas/lapack
**               connector.
**
**	 Author: Brian Amberg, 2007
**		 Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <iomanip>
#include "tMatrix.h"
#include "Exception.h"
#include "t2Matrix.h"
#include "t3Matrix.h"
#include "t4Matrix.h"


namespace gravis
{
  ////////////////////////////////////////////////////////////////////
  // The Matrix classes. We distinguish between vectors and matrices
  // and allow for views of full matrices. No spacing magic is used, these are
  // simple dense matrices in column-first order

  template <class T> class tVarVector;
  template <class T> class tVectorView;
  template <class T> class tConstVectorView;

  template <class T> class tVarMatrix;
  template <class T> class tMatrixView;
  template <class T> class tConstMatrixView;

  //
  // IMPLEMENTATION
  //

  namespace matrix
  {
    /**
     * Set all to zero
     **/
    template <class VectorOrMatrix>
    inline static
    void clear(VectorOrMatrix& v)
    {
      if (v.size()>0) memset( &v[0], 0, sizeof(v[0])*v.size());
    }

    /**
     * Fill with equal elements
     **/
    template <class VectorOrMatrix>
    inline static
    void fill(VectorOrMatrix& v, const typename VectorOrMatrix::scalar& value)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] = value;
    }

    /**
     * Clamp all elements of a matrix
     **/
    template <class VectorOrMatrix> static inline void clamp(VectorOrMatrix& v, const typename VectorOrMatrix::scalar& min, const typename VectorOrMatrix::scalar& max)
    {
      for (size_t i=0; i<v.size(); ++i)
      {
        if (v[i]<min) v[i] = min;
        if (max<v[i]) v[i] = max;
      }
    }

    namespace priv
    {

      //////////////////////// MEMORY HELPER, HAD A HOOK FOR TRACKING ALLOCATIONS, BUT REMOVED FOR LIBGRAVIS //////
      template <class T>
      inline static
      T* alloc_arr(const std::string& title, const size_t h, const size_t w=1)
      {
        T* r = new T[h*w];
        return r;
      }

      template <class T>
      inline static
      void free_arr(const T* p)
      {
        delete [] p;
      }

      template <class T>
      inline static
      void copy_arr(T* t, const T* s, const size_t sz)
      {
        if (sz>0)
          memcpy(t, s, sz*sizeof(T));
      }

      /////////////////////// HELPER //////////////////////

      /**
       * Clamp a value
       **/
      template <class T>
      static inline
      void clamp(T& v, const T& min, const T& max)
      {
        if (v<min) v = min;
        if (max<v) v = max;
      }

      /**
       * Take the square of a value
       **/
      template <class T> inline T sqr(const T& v)
      {
        return v*v;
      }
    }
  }


  /**
   * A thin c++ matrix wrapper around a slice of memory
   *
   * These matrix classes allow easy access of the blas/lapack functions.
   *
   * They are relatively rough, as they try to be as simple as possible. In my
   * view it is not good to make c++ behave like matlab, as the only advantage
   * of c++ over matlab is more control. These classes give maximum control.
   **/
  template <class T>
  class tVectorView
  {
    public:
      typedef T scalar;

      size_t h;
      T* const data;

      tVectorView(T* data, size_t h) : h(h), data(data)  {    }
      /**
       * Create a view of the other matrix
       **/
      tVectorView(tVectorView<T> &o) : h(o.h), data(o.data) {}
      tVectorView(tVarVector<T>  &o) : h(o.h), data(o.data) {}
      template <size_t mh>
      tVectorView(tMatrix<T, mh, 1> &m) : h(mh),  data(&m[0]) {}
      /**
       * Copy another vector into this vector.
       **/
      tVectorView& operator=(const tConstVectorView<T> &o)
      {
        GRAVIS_CHECK(o.h==h, "Incompatible size");
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
        return *this;
      }
      tVectorView& operator=(const tVectorView<T>      &o)
      {
        GRAVIS_CHECK(o.h==h, "Incompatible size");
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
        return *this;
      }
      tVectorView& operator=(const tVarVector<T>       &o)
      {
        GRAVIS_CHECK(o.h==h, "Incompatible size");
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
        return *this;
      }

      inline       T& operator[](size_t i)
      {
        checkAccess1( i, h );
        return data[i];
      };
      inline const T& operator[](size_t i) const
      {
        checkAccess1( i, h );
        return data[i];
      };

      inline       T& operator()(size_t i)
      {
        checkAccess1( i, h );
        return data[i];
      };
      inline const T& operator()(size_t i) const
      {
        checkAccess1( i, h );
        return data[i];
      };

      inline       T& clampedAccess(int i)
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        return operator()(i);
      }
      inline const T& clampedAccess(int i) const
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        return operator()(i);
      }

      inline size_t size() const
      {
        return h;
      }

      /**
       * Convenience function to clear a matrix
       **/
      inline void clear()
      {
        gravis::matrix::clear(*this);
      }
      /**
       * Convenience functions to fill a matrix
       **/
      inline void fill(const T& e)
      {
        gravis::matrix::fill(*this, e);
      }
      /**
       * Convenience functions to clamp all elements of a matrix
       **/
      inline void clamp(const T& min, const T& max)
      {
        gravis::matrix::clamp(*this, min, max);
      }
  };

  /**
   * A thin c++ matrix wrapper around a slice of memory
   *
   * These matrix classes allow easy access of the blas/lapack functions.
   *
   * They are relatively rough, as they try to be as simple as possible. In my
   * view it is not good to make c++ behave like matlab, as the only advantage
   * of c++ over matlab is more control. These classes give maximum control.
   **/
  template <class T>
  class tConstVectorView
  {
    public:
      typedef T scalar;
      size_t h;
      const T* const data;

      tConstVectorView(const T* data, size_t h) : h(h),   data(data) {}
      tConstVectorView(const tConstVectorView& o) : h(o.h), data(o.data) {}
      tConstVectorView(const tVectorView<T>   &o) : h(o.h), data(o.data) {}
      tConstVectorView(const tVarVector<T>    &o) : h(o.h), data(o.data) {}
      template <size_t mh>
      tConstVectorView(const tMatrix<T, mh, 1> &m) : h(mh), data(&m[0])  {}

      inline const T& operator[](size_t i) const
      {
        checkAccess1( i, h );
        return data[i];
      }
      inline const T& operator()(size_t i) const
      {
        checkAccess1( i, h );
        return data[i];
      }

      inline const T& clampedAccess(int i) const
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        return operator()(i);
      }

      inline size_t size() const
      {
        return h;
      }
  };

  /**
   * A matrix with memory allocated on the heap.
   *
   * The semantic of operations on this vector is different from the vector
   * views. Assigning something to this vector is a copy operation, while for the
   * views it is just a pointer assignment.
   *
   * They are relatively rough, as they try to be as simple as possible. In my
   * view it is not good to make c++ behave like matlab, as the only advantage
   * of c++ over matlab is more control. These classes give maximum control.
   **/
  template <class T>
  class tVarVector
  {
    public:
      typedef T scalar;

      size_t h;
      T* data;

      std::string title;

      tVarVector(size_t h, const std::string& title="UNNAMED:VECTOR") : h(h), data(matrix::priv::alloc_arr<T>(title, h)), title(title) {}
      tVarVector(const std::string& title="UNNAMED:VECTOR") : h(0), data(matrix::priv::alloc_arr<T>(title, h)), title(title) {}
      ~tVarVector()
      {
        matrix::priv::free_arr(data);
      };
      /**
       * Copy another vector into this vector
       **/
      tVarVector(const tConstVectorView<T> &o, const std::string& title="UNNAMED:VECTOR") : h(o.h), data(matrix::priv::alloc_arr<T>(title, h)), title(title)
      {
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
      }
      tVarVector(const tVectorView<T>      &o, const std::string& title="UNNAMED:VECTOR") : h(o.h), data(matrix::priv::alloc_arr<T>(title, h)), title(title)
      {
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
      }
      tVarVector(const tVarVector&          o, const std::string& title="UNNAMED:VECTOR") : h(o.h), data(matrix::priv::alloc_arr<T>(title, h)), title(title)
      {
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
      }
      template <size_t mh>
      tVarVector(const tMatrix<T, mh, 1>   &o, const std::string& title="UNNAMED:VECTOR") : h(mh), data(matrix::priv::alloc_arr<T>(title, h)), title(title)
      {
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
      }
      /**
       * Copy another vector into this vector.
       * Will loose old data reference, beware.
       **/
      tVarVector& operator=(const tConstVectorView<T> &o)
      {
        resize(o.h);
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
        return *this;
      }
      tVarVector& operator=(const tVectorView<T>      &o)
      {
        resize(o.h);
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
        return *this;
      }
      tVarVector& operator=(const tVarVector<T>       &o)
      {
        resize(o.h);
        //for (size_t i=0; i<h; ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, h);
        return *this;
      }

      /**
       * Will loose old data reference, beware.
       **/
      void resize(size_t h)
      {
        if (h == this->h)
          return;
        if (h > this->h)
        {
          T* new_data = matrix::priv::alloc_arr<T>(title, h);
          std::swap(data, new_data);
          matrix::priv::free_arr(new_data);
        }
        this->h = h;
      }

      inline       T& operator[](size_t i)
      {
        checkAccess1( i, h );
        return data[i];
      };
      inline const T& operator[](size_t i) const
      {
        checkAccess1( i, h );
        return data[i];
      };

      inline       T& operator()(size_t i)
      {
        checkAccess1( i, h );
        return data[i];
      };
      inline const T& operator()(size_t i) const
      {
        checkAccess1( i, h );
        return data[i];
      };

      inline       T& clampedAccess(int i)
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        return operator()(i);
      }
      inline const T& clampedAccess(int i) const
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        return operator()(i);
      }

      inline size_t size() const
      {
        return h;
      }

      /**
       * Convenience function to clear a matrix
       **/
      inline void clear()
      {
        gravis::matrix::clear(*this);
      }
      /**
       * Convenience functions to fill a matrix
       **/
      inline void fill(const T& e)
      {
        gravis::matrix::fill(*this, e);
      }
      /**
       * Convenience functions to clamp all elements of a matrix
       **/
      inline void clamp(const T& min, const T& max)
      {
        gravis::matrix::clamp(*this, min, max);
      }
  };

  /**
   * A thin c++ matrix wrapper around a slice of memory
   **/
  template <class T>
  class tMatrixView
  {
    public:
      typedef T scalar;
      size_t h, w;
      T* const data;

      tMatrixView(T* data, size_t h, size_t w) : h(h), w(w), data(data) {}
      tMatrixView(tMatrixView<T> &o) : h(o.h), w(o.w), data(o.data) {}
      tMatrixView(tVarMatrix<T> &o)  : h(o.h), w(o.w), data(o.data) {}
      template <size_t mw, size_t mh>
      tMatrixView(tMatrix<T, mw, mh> &m) : h(mh), w(mw), data(&m[0]) {}
      /**
       * Copy another vector into this vector.
       **/
      tMatrixView& operator=(const tConstMatrixView<T> &o)
      {
        GRAVIS_CHECK(o.h==h && o.w==w, "Incompatible size");
        //for (size_t i=0; i<size(); ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, size());
        return *this;
      }
      tMatrixView& operator=(const tMatrixView<T>      &o)
      {
        GRAVIS_CHECK(o.h==h && o.w==w, "Incompatible size");
        //for (size_t i=0; i<size(); ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, size());
        return *this;
      }
      tMatrixView& operator=(const tVarMatrix<T>       &o)
      {
        GRAVIS_CHECK(o.h==h && o.w==w, "Incompatible size");
        //for (size_t i=0; i<size(); ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, size());
        return *this;
      }

      inline       T& operator[](size_t i)
      {
        checkAccess1( i, h*w );
        return data[i];
      };
      inline const T& operator[](size_t i) const
      {
        checkAccess1( i, h*w );
        return data[i];
      };

      inline       T& operator()(size_t i, size_t j)
      {
        checkAccess2( i,j, h, w );
        return data[i + j*h];
      };
      inline const T& operator()(size_t i, size_t j) const
      {
        checkAccess2( i,j, h, w );
        return data[i + j*h];
      };

      inline       T& clampedAccess(int i, int j)
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        matrix::priv::clamp(j, 0, int(w)-1);
        return operator()(i,j);
      }
      inline const T& clampedAccess(int i, int j) const
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        matrix::priv::clamp(j, 0, int(w)-1);
        return operator()(i,j);
      }

      inline size_t size() const
      {
        return h*w;
      }

      /**
       * Convenience function to clear a matrix
       **/
      inline void clear()
      {
        gravis::matrix::clear(*this);
      }
      /**
       * Convenience functions to fill a matrix
       **/
      inline void fill(const T& e)
      {
        gravis::matrix::fill(*this, e);
      }
      /**
       * Convenience functions to clamp all elements of a matrix
       **/
      inline void clamp(const T& min, const T& max)
      {
        gravis::matrix::clamp(*this, min, max);
      }
  };

  /**
   * A thin c++ matrix wrapper around a slice of memory
   **/
  template <class T>
  class tConstMatrixView
  {
    public:
      typedef T scalar;
      size_t h, w;
      const T* const data;

      tConstMatrixView(const T* data, size_t h, size_t w) : h(h), w(w), data(data) {}
      tConstMatrixView(const tConstMatrixView& o) : h(o.h), w(o.w), data(o.data) {}
      tConstMatrixView(const tMatrixView<T> &o)   : h(o.h), w(o.w), data(o.data) {}
      tConstMatrixView(const tVarMatrix<T> &o)    : h(o.h), w(o.w), data(o.data) {}
      template <size_t mw, size_t mh>
      tConstMatrixView(const tMatrix<T, mw, mh> &m) : h(mh), w(mw), data(&m[0]) {}

      inline const T& operator[](size_t i          ) const
      {
        checkAccess1( i, h*w );
        return data[i];
      }
      inline const T& operator()(size_t i, size_t j) const
      {
        checkAccess2( i,j, h, w );
        return data[i + j*h];
      }

      inline const T& clampedAccess(int i, int j) const
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        matrix::priv::clamp(j, 0, int(w)-1);
        return operator()(i,j);
      }

      inline size_t size() const
      {
        return h*w;
      }
  };

  /**
   * A matrix with memory allocated on the heap
   **/
  template <class T>
  class tVarMatrix
  {
    public:
      typedef T scalar;

      size_t h,w;
      T* data;

      std::string title;

      tVarMatrix(size_t h, size_t w, const std::string& title="UNNAMED:MATRIX") : h(h), w(w), data(matrix::priv::alloc_arr<T>(title, h, w)), title(title) {}
      tVarMatrix(const std::string& title="UNNAMED:MATRIX") : h(0), w(0), data(matrix::priv::alloc_arr<T>(title, h, w)), title(title) {}
      ~tVarMatrix()
      {
        matrix::priv::free_arr(data);
      };
      /**
       * Copy another matrix into this matrix
       **/
      tVarMatrix(const tConstMatrixView<T> &o, const std::string& title="UNNAMED:MATRIX") : h(o.h), w(o.w), data(matrix::priv::alloc_arr<T>(title, o.h, o.w)), title(title)
      {
        matrix::priv::copy_arr(data, o.data, size());
      }
      /**
       * Copy another matrix into this matrix
       **/
      tVarMatrix(const tMatrixView<T>      &o, const std::string& title="UNNAMED:MATRIX") : h(o.h), w(o.w), data(matrix::priv::alloc_arr<T>(title, o.h, o.w)), title(title)
      {
        matrix::priv::copy_arr(data, o.data, size());
      }
      /**
       * Copy another matrix into this matrix
       **/
      tVarMatrix(const tVarMatrix<T>       &o, const std::string& title="UNNAMED:MATRIX") : h(o.h), w(o.w), data(matrix::priv::alloc_arr<T>(title, o.h, o.w)), title(title)
      {
        matrix::priv::copy_arr(data, o.data, size());
      }
      /**
       * Copy another matrix into this matrix
       **/
      template <size_t mh, size_t mw>
      tVarMatrix(const tMatrix<T, mh, mw>  &o, const std::string& title="UNNAMED:MATRIX") : h(mh), w(mw), data(matrix::priv::alloc_arr<T>(title, o.h, o.w)), title(title)
      {
        matrix::priv::copy_arr(data, o.data, size());
      }
      /**
       * Copy another matrix into this matrix
       **/
      tVarMatrix(const t2Matrix<T>  &o, const std::string& title="UNNAMED:MATRIX") : h(2), w(2), data(matrix::priv::alloc_arr<T>(title, h, w)), title(title)
      {
        matrix::priv::copy_arr(data, o.m, size());
      }
      /**
       * Copy another matrix into this matrix
       **/
      tVarMatrix(const t3Matrix<T>  &o, const std::string& title="UNNAMED:MATRIX") : h(3), w(3), data(matrix::priv::alloc_arr<T>(title, h, w)), title(title)
      {
        matrix::priv::copy_arr(data, o.m, size());
      }
      /**
       * Copy another matrix into this matrix
       **/
      tVarMatrix(const t4Matrix<T>  &o, const std::string& title="UNNAMED:MATRIX") : h(4), w(4), data(matrix::priv::alloc_arr<T>(title, h, w)), title(title)
      {
        matrix::priv::copy_arr(data, o.m, size());
      }
      /**
       * Copy another matrix into this matrix.
       * Will loose old data reference, beware.
       **/
      tVarMatrix& operator=(const tConstMatrixView<T> &o)
      {
        resize(o.h,o.w);
        //for (size_t i=0; i<size(); ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, size());
        return *this;
      }
      tVarMatrix& operator=(const tMatrixView<T>      &o)
      {
        resize(o.h,o.w);
        //for (size_t i=0; i<size(); ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, size());
        return *this;
      }
      tVarMatrix& operator=(const tVarMatrix<T>       &o)
      {
        resize(o.h,o.w);
        //for (size_t i=0; i<size(); ++i) (*this)[i] = o[i];
        matrix::priv::copy_arr(data, o.data, size());
        return *this;
      }

      /**
       * Will loose old data reference, beware.
       **/
      void resize(size_t h, size_t w)
      {
        if (h == this->h && w == this->w)
          return;
        if (h*w>size())
        {
          T* new_data = matrix::priv::alloc_arr<T>(title, h, w);
          std::swap(data, new_data);
          matrix::priv::free_arr(new_data);
        }
        this->h = h;
        this->w = w;
      }

      inline       T& operator[](size_t i)
      {
        checkAccess1( i, h*w );
        return data[i];
      };
      inline const T& operator[](size_t i) const
      {
        checkAccess1( i, h*w );
        return data[i];
      };

      inline       T& operator()(size_t i, size_t j)
      {
        checkAccess2( i, j, h, w );
        return data[i+j*h];
      };
      inline const T& operator()(size_t i, size_t j) const
      {
        checkAccess2( i, j, h, w );
        return data[i+j*h];
      };

      inline       T& clampedAccess(int i, int j)
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        matrix::priv::clamp(j, 0, int(w)-1);
        return operator()(i,j);
      }
      inline const T& clampedAccess(int i, int j) const
      {
        matrix::priv::clamp(i, 0, int(h)-1);
        matrix::priv::clamp(j, 0, int(w)-1);
        return operator()(i,j);
      }

      inline size_t size() const
      {
        return h*w;
      }

      /**
       * Convenience function to clear a matrix
       **/
      inline void clear()
      {
        gravis::matrix::clear(*this);
      }
      /**
       * Convenience functions to fill a matrix
       **/
      inline void fill(const T& e)
      {
        gravis::matrix::fill(*this, e);
      }
      /**
       * Convenience functions to clamp all elements of a matrix
       **/
      inline void clamp(const T& min, const T& max)
      {
        gravis::matrix::clamp(*this, min, max);
      }
  };

  /**
   * Matrix and vector operations
   **/
  namespace matrix
  {

    template <class T>
    inline static
    void display( const tConstVectorView<T> &v)
    {
      std::cout << "Vector: " << v.h << std::endl;
      for (size_t i=0; i<v.h; ++i)
      {
        std::cout << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(3) << v[i] << std::endl;
      }
    }

    template <class T>
    inline static
    void display( const tVarVector<T> &v)
    {
      display<T>( tConstVectorView<T>(v) );
    }

    template <class T>
    inline static
    void display( const tVectorView<T> &v)
    {
      display<T>( tConstVectorView<T>(v) );
    }

    template <class T>
    inline static
    void display( const tConstMatrixView<T> &v)
    {
      std::cout << "Matrix: " << v.h << "x" << v.w << std::endl;
      for (size_t i=0; i<v.h; ++i)
      {
        for (size_t j=0; j<v.w; ++j)
        {
          std::cout << std::fixed << std::right << std::setw(8) << std::showpoint << std::setprecision(3) << v(i,j) << " ";
        }
        std::cout << std::endl;
      }
    }

    template <class T>
    inline static
    void display( const tVarMatrix<T> &v)
    {
      display<T>( tConstMatrixView<T>(v) );
    }

    template <class T>
    inline static
    void display( const tMatrixView<T> &v)
    {
      display<T>( tConstMatrixView<T>(v) );
    }

    /**
     * Find the largest element
     **/
    template <class VectorOrMatrix>
    inline static
    typename VectorOrMatrix::scalar max(const VectorOrMatrix& v)
    {
      size_t mi = 0;
      for (size_t i=1; i<v.size(); ++i) if (v[i] > v[mi]) mi = i;
      return v[mi];
    }

    /**
     * Find the smallest element
     **/
    template <class VectorOrMatrix>
    inline static
    typename VectorOrMatrix::scalar min(const VectorOrMatrix& v)
    {
      size_t mi = 0;
      for (size_t i=1; i<v.size(); ++i) if (v[i] < v[mi]) mi = i;
      return v[mi];
    }

    /** Arithmethic operations with scalars -= **/
    template <class T>                     inline static void sub(tVectorView<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] -= s;
    }
    /** Arithmethic operations with scalars -= **/
    template <class T>                     inline static void sub(tVarVector<T>   &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] -= s;
    }
    /** Arithmethic operations with scalars -= **/
    template <class T>                     inline static void sub(tMatrixView<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] -= s;
    }
    /** Arithmethic operations with scalars -= **/
    template <class T>                     inline static void sub(tVarMatrix<T>   &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] -= s;
    }
    /** Arithmethic operations with scalars -= **/
    template <class T, size_t h, size_t w> inline static void sub(tMatrix<T,h,w>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] -= s;
    }

    /** Arithmethic operations with scalars += **/
    template <class T>                     inline static void add(tVectorView<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] += s;
    }
    /** Arithmethic operations with scalars += **/
    template <class T>                     inline static void add(tVarVector<T>   &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] += s;
    }
    /** Arithmethic operations with scalars += **/
    template <class T>                     inline static void add(tMatrixView<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] += s;
    }
    /** Arithmethic operations with scalars += **/
    template <class T>                     inline static void add(tVarMatrix<T>   &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] += s;
    }
    /** Arithmethic operations with scalars += **/
    template <class T, size_t h, size_t w> inline static void add(tMatrix<T,h,w>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] += s;
    }

    /** Arithmethic operations with scalars *= **/
    template <class T>                     inline static void mult(tVectorView<T> &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] *= s;
    }
    /** Arithmethic operations with scalars *= **/
    template <class T>                     inline static void mult(tVarVector<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] *= s;
    }
    /** Arithmethic operations with scalars *= **/
    template <class T>                     inline static void mult(tMatrixView<T> &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] *= s;
    }
    /** Arithmethic operations with scalars *= **/
    template <class T>                     inline static void mult(tVarMatrix<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] *= s;
    }
    /** Arithmethic operations with scalars *= **/
    template <class T, size_t h, size_t w> inline static void mult(tMatrix<T,h,w> &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] *= s;
    }

    /** Arithmethic operations with scalars /= **/
    template <class T>                     inline static void div(tVectorView<T> &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] /= s;
    }
    /** Arithmethic operations with scalars /= **/
    template <class T>                     inline static void div(tVarVector<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] /= s;
    }
    /** Arithmethic operations with scalars /= **/
    template <class T>                     inline static void div(tMatrixView<T> &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] /= s;
    }
    /** Arithmethic operations with scalars /= **/
    template <class T>                     inline static void div(tVarMatrix<T>  &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] /= s;
    }
    /** Arithmethic operations with scalars /= **/
    template <class T, size_t h, size_t w> inline static void div(tMatrix<T,h,w> &v, const T& s)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] /= s;
    }

    /** Arithmethic operations with scalars A=-A **/
    template <class T>                     inline static void negate(tVectorView<T>  &v)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] = -v[i];
    }
    /** Arithmethic operations with scalars A=-A **/
    template <class T>                     inline static void negate(tVarVector<T>  &v)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] = -v[i];
    }
    /** Arithmethic operations with scalars A=-A **/
    template <class T>                     inline static void negate(tMatrixView<T>  &v)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] = -v[i];
    }
    /** Arithmethic operations with scalars A=-A **/
    template <class T>                     inline static void negate(tVarMatrix<T>  &v)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] = -v[i];
    }


    /** Arithmethic operations with matrices. Element wise addition += **/
    template <class T>                     inline static void add(tVectorView<T>  &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] += v2[i];
    }
    /** Arithmethic operations with matrices. Element wise addition += **/
    template <class T>                     inline static void add(tVarVector<T>   &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] += v2[i];
    }
    /** Arithmethic operations with matrices. Element wise addition += **/
    template <class T>                     inline static void add(tMatrixView<T>  &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] += v2[i];
    }
    /** Arithmethic operations with matrices. Element wise addition += **/
    template <class T>                     inline static void add(tVarMatrix<T>   &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] += v2[i];
    }
    /** Arithmethic operations with matrices. Element wise addition += **/
    template <class T, size_t h, size_t w> inline static void add(tMatrix<T,h,w>  &v, const tMatrix<T,h,w>&      v2)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] += v2[i];
    }


    /** Arithmethic operations with matrices. Element wise subtraction -= **/
    template <class T>                     inline static void sub(tVectorView<T>  &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] -= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise subtraction -= **/
    template <class T>                     inline static void sub(tVarVector<T>   &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] -= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise subtraction -= **/
    template <class T>                     inline static void sub(tMatrixView<T>  &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size, "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] -= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise subtraction -= **/
    template <class T>                     inline static void sub(tVarMatrix<T>   &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] -= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise subtraction -= **/
    template <class T, size_t h, size_t w> inline static void sub(tMatrix<T,h,w>  &v, const tMatrix<T,h,w>&      v2)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] -= v2[i];
    }

    /** Arithmethic operations with matrices. Element wise multiplication *= **/
    template <class T>                     inline static void elmul(tVectorView<T>  &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] *= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise multiplication *= **/
    template <class T>                     inline static void elmul(tVarVector<T>   &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] *= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise multiplication *= **/
    template <class T>                     inline static void elmul(tMatrixView<T>  &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] *= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise multiplication *= **/
    template <class T>                     inline static void elmul(tVarMatrix<T>   &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() != v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] *= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise multiplication *= **/
    template <class T, size_t h, size_t w> inline static void elmul(tMatrix<T,h,w>  &v, const tMatrix<T,h,w>&      v2)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] *= v2[i];
    }

    /** Arithmethic operations with matrices. Element wise division /= **/
    template <class T>                     inline static void eldiv(tVectorView<T>  &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() == v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] /= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise division /= **/
    template <class T>                     inline static void eldiv(tVarVector<T>   &v, const tConstVectorView<T>& v2)
    {
      GRAVIS_CHECK(v.size() == v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] /= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise division /= **/
    template <class T>                     inline static void eldiv(tMatrixView<T>  &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() == v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] /= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise division /= **/
    template <class T>                     inline static void eldiv(tVarMatrix<T>   &v, const tConstMatrixView<T>& v2)
    {
      GRAVIS_CHECK(v.size() == v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] /= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise division /= **/
    template <class T>                     inline static void eldiv(tVarMatrix<T>   &v, const tVarMatrix<T>&       v2)
    {
      GRAVIS_CHECK(v.size() == v2.size(), "Size must be equal for addition");
      for(size_t i=0; i<v.size(); ++i) v[i] /= v2[i];
    }
    /** Arithmethic operations with matrices. Element wise division /= **/
    template <class T, size_t h, size_t w> inline static void eldiv(tMatrix<T,h,w>  &v, const tMatrix<T,h,w>&      v2)
    {
      for(size_t i=0; i<v.size(); ++i) v[i] /= v2[i];
    }

    /**
     * Per element larger than test for scalars
     **/
    template <class OutMatrix, class InMatrix>
    inline static
    void cmpLarger( OutMatrix& Mout, const InMatrix& M, const typename InMatrix::scalar& t)
    {
      GRAVIS_CHECK(Mout.h == M.h && Mout.w == M.w, "Incompatible sizes");
      const int S=Mout.size();
      for (int i=0; i<S; ++i) Mout[i] = M[i]>t ? '\xFF' : '\x00';
    }

    /**
     * Per element smaller than test for scalars
     **/
    template <class OutMatrix, class InMatrix>
    inline static
    void cmpSmaller( OutMatrix& Mout, const InMatrix& M, const typename InMatrix::scalar& t)
    {
      GRAVIS_CHECK(Mout.h == M.h && Mout.w == M.w, "Incompatible sizes");
      int i;
      const int S=Mout.size();
      for (int i=0; i<S; ++i) Mout[i] = M[i]<t ? '\xFF' : '\x00';
    }

    /**
     * Per element equality test for scalars
     **/
    template <class OutMatrix, class InMatrix>
    inline static
    void cmpEqual( OutMatrix& Mout, const InMatrix& M, const typename InMatrix::scalar& t)
    {
      GRAVIS_CHECK(Mout.h == M.h && Mout.w == M.w, "Incompatible sizes");
      int i;
      const int S=Mout.size();
      for (int i=0; i<S; ++i) Mout[i] = M[i]==t ? '\xFF' : '\x00';
    }

    /**
     * Inset one matrix into another
     **/
    template <class OutMatrix, class InMatrix>
    inline static
    void inset( OutMatrix& Out, const InMatrix& In, const size_t row, const size_t col=0)
    {
      if ((In.h == 0) || (In.w == 0) || (col>=Out.w) || (row>=Out.h)) return;
      size_t h=std::min(In.h, Out.h-row);
      size_t w=std::min(In.w, Out.w-col);
      for (size_t j=0; j<w; ++j)
      {
        memcpy( &Out(row, col+j), &In(0, j), sizeof(In[0])*h );
      }
    }

    /**
     * Inset one matrix into another
     **/
    template <class T>
    inline static
    void inset( tVarVector<T> &Out, const tConstVectorView<T> &In, const size_t row)
    {
      if ((In.h == 0) || (row>=Out.h)) return;
      size_t h=std::min(In.h, Out.h-row);
      memcpy( &Out(row, 0), &In(0, 0), sizeof(In[0])*h );
    }

    /**
     * Inset one matrix into another
     **/
    template <class T>
    inline static
    void inset( tVarVector<T> &Out, const tVarVector<T> &In, const size_t row)
    {
      if ((In.h == 0) || (row>=Out.h)) return;
      size_t h=std::min(In.h, Out.h-row);
      memcpy( &Out[row], &In[0], sizeof(In[0])*h );
    }

    /**
     * Matrix Convolution
     *
     * TODO: Do not use checked access in the main region of the image, use it
     * only on the borders
     **/
    template <class OutMatrix, class InMatrixImg, class InMatrixMask>
    inline static
    void conv2( OutMatrix& Iout, const InMatrixImg& I, const InMatrixMask& F )
    {
      GRAVIS_CHECK( Iout.w == I.w && Iout.h == I.h, "Matrix sizes are not compatible" );
      Iout.clear();
      const int ox(F.w/2);
      const int oy(F.h/2);
      const int W=I.w;
      int j;
#ifdef _OPENMP
      #pragma omp parallel for default(none) private(j) shared(I,Iout,F)
#endif
      for (j=0; j<W; ++j)
        for (int i=0; i<int(I.h); ++i)
          for (int l=0; l<int(F.w); ++l)
            for (int k=0; k<int(F.h); ++k)
              Iout(i,j) += F(k,l) * I.clampedAccess(i+oy-k,j+ox-l);
    }

    /**
     * Erosion of a binary matrix
     *
     * TODO: The outmost border is not handled
     **/
    template <class Matrix>
    inline static
    void erode( Matrix& m)
    {
      const tVarMatrix<typename Matrix::scalar> M(m); // Make a copy
      int j;
      const int W=M.w;
#ifdef _OPENMP
      #pragma omp parallel for default(none) private(j) shared(m)
#endif
      for (j=1; j<W-1; ++j)
        for (size_t i=1; i<M.h-1; ++i)
          m(i,j) = m(i,j) && M(i-1,j) && M(i,j-1) && M(i+1,j) && M(i,j+1);
    }


    template <class OutMatrix, class InMatrix1, class InMatrix2>
    inline static
    void mult_elementwise(OutMatrix& m, const InMatrix1& m1, const InMatrix2& m2)
    {
      GRAVIS_CHECK( m.size() == m1.size(), "Matrix sizes incompatible");
      GRAVIS_CHECK( m.size() == m2.size(), "Matrix sizes incompatible");
      const size_t s = m.size();
      for (size_t i=0; i<s; ++i) m[i] = m1[i] * m2[i];
    }

    /** Sum up the values of a matrix *= **/
    template <class T> inline static T sum(tConstMatrixView<T> &v)
    {
      if (v.size()==0) return T(0);
      T r=v[0];
      for(size_t i=1; i<v.size(); ++i) r += v[i];
      return r;
    }
    /** Sum up the values of a matrix *= **/
    template <class T> inline static T sum(tConstVectorView<T> &v)
    {
      if (v.size()==0) return T(0);
      T r=v[0];
      for(size_t i=1; i<v.size(); ++i) r += v[i];
      return r;
    }
    /** Sum up the values of a matrix *= **/
    template <class T> inline static T sum(tMatrixView<T> &v)
    {
      if (v.size()==0) return T(0);
      T r=v[0];
      for(size_t i=1; i<v.size(); ++i) r += v[i];
      return r;
    }
    /** Sum up the values of a matrix *= **/
    template <class T> inline static T sum(tVectorView<T> &v)
    {
      if (v.size()==0) return T(0);
      T r=v[0];
      for(size_t i=1; i<v.size(); ++i) r += v[i];
      return r;
    }
    /** Sum up the values of a matrix *= **/
    template <class T> inline static T sum(tVarMatrix<T> &v)
    {
      if (v.size()==0) return T(0);
      T r=v[0];
      for(size_t i=1; i<v.size(); ++i) r += v[i];
      return r;
    }
    /** Sum up the values of a matrix *= **/
    template <class T> inline static T sum(tVarVector<T> &v)
    {
      if (v.size()==0) return T(0);
      T r=v[0];
      for(size_t i=1; i<v.size(); ++i) r += v[i];
      return r;
    }

    // Matrix input output
    template <class T>
    static
    void load(tVarMatrix<T> &M, const std::string& fn)
    {
      char mmid0[33] = "GRAVIS_VAR_MATRIX               ";
      char mmid1[33] = "GRAVIS_VAR_MATRIX               ";
      std::ifstream stream(fn.c_str(), std::ifstream::binary);
      uint8_t uint32_size;
      uint8_t T_size;
      uint32_t h,w;
      uint16_t endianness;
      stream.read(mmid1, 32);
      stream.read((char*)&endianness, 2);
      stream.read((char*)&uint32_size, 1);
      stream.read((char*)&T_size, 1);
      stream.read((char*)&h, sizeof(h));
      stream.read((char*)&w, sizeof(w));
      GRAVIS_CHECK( 0==strncmp( mmid0, mmid1, 31 ),"Not a gravis var matrix file" );
      GRAVIS_CHECK( endianness  == 0x0001,         "Wrong endianness");
      GRAVIS_CHECK( uint32_size == 4,              "Wrong size_t size");
      GRAVIS_CHECK( T_size      == sizeof(T),      "Wrong type in matrix file");
      M.resize(h,w);
      stream.read((char*)M.data, sizeof(T)*M.size());
    }

    template <class T>
    static
    void save(const std::string& fn, const tConstMatrixView<T> &v)
    {
      char mmid[33] = "GRAVIS_VAR_MATRIX               ";
      std::ofstream stream(fn.c_str(), std::ofstream::binary);
      uint8_t uint32_size = sizeof(uint32_t);
      uint8_t T_size      = sizeof(T);
      uint32_t h = v.h, w = v.w;
      uint16_t endianness = 0x0001;
      stream.write(mmid, 32);
      stream.write((char*)&endianness, 2);
      stream.write((char*)&uint32_size, 1);
      stream.write((char*)&T_size, 1);
      stream.write((char*)&h, sizeof(h));
      stream.write((char*)&w, sizeof(w));
      stream.write((char*)v.data, sizeof(T)*v.size());
    }

    template <class T>
    static
    void load(tVarVector<T> &v, const std::string& fn)
    {
      char mmid0[33] = "GRAVIS_VAR_VECTOR              ";
      char mmid1[33] = "GRAVIS_VAR_VECTOR              ";
      std::ifstream stream(fn.c_str(), std::ifstream::binary);
      uint8_t uint32_size;
      uint8_t T_size;
      uint32_t k;
      uint16_t endianness;
      stream.read(mmid1, 32);
      stream.read((char*)&endianness, 2);
      stream.read((char*)&uint32_size, 1);
      stream.read((char*)&T_size, 1);
      stream.read((char*)&k, sizeof(k));
      GRAVIS_CHECK( 0 == strncmp( mmid0, mmid1, 31 ), "Not a gravis var vector file" );
      GRAVIS_CHECK( endianness  == 0x0001,            "Wrong endianness");
      GRAVIS_CHECK( uint32_size == 4,                 "Wrong uint32 size");
      GRAVIS_CHECK( T_size      == sizeof(T),         "Wrong type in model file");
      v.resize(k);
      stream.read((char*)v.data, sizeof(T)*v.size());
    }

    template <class T>
    static
    void save(const std::string& fn, const tConstVectorView<T> &v)
    {
      char mmid[33] = "GRAVIS_VAR_VECTOR               ";
      std::ofstream stream(fn.c_str(), std::ofstream::binary);
      uint8_t uint32_size = sizeof(uint32_t);
      uint8_t T_size      = sizeof(T);
      uint16_t endianness = 0x0001;
      uint32_t k = v.size();
      stream.write(mmid, 32);
      stream.write((char*)&endianness, 2);
      stream.write((char*)&uint32_size, 1);
      stream.write((char*)&T_size, 1);
      stream.write((char*)&k, sizeof(k));
      stream.write((char*)v.data, sizeof(T)*v.size());
    }

    template <class T> static inline void clamp(tMatrixView<T> &v, const T& min, const T& max)
    {
      for (size_t i=0; i<v.size(); ++i) priv::clamp(v[i], min, max);
    }
    template <class T> static inline void clamp(tVarMatrix <T> &v, const T& min, const T& max)
    {
      for (size_t i=0; i<v.size(); ++i) priv::clamp(v[i], min, max);
    }
    template <class T> static inline void clamp(tVectorView<T> &v, const T& min, const T& max)
    {
      for (size_t i=0; i<v.size(); ++i) priv::clamp(v[i], min, max);
    }
    template <class T> static inline void clamp(tVarVector <T> &v, const T& min, const T& max)
    {
      for (size_t i=0; i<v.size(); ++i) priv::clamp(v[i], min, max);
    }

  }

  /**
   * Read Variable size matrices from a stream
   **/
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, tVectorView<T>& arg)
  {
    size_t h = arg.h;
    std::string t;
    is >> t;
    if (t != "[") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A vector should start with [", t);
    for (size_t j=0; j<h; ++j)
      is >> arg[j];
    is >> t;
    if (t != "]") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A vector should end with ]", t);
    return is;
  }
  /**
   * Read Variable size matrices from a stream
   **/
  template <class T>
  inline
  std::istream& operator>> (std::istream& is, tVarVector<T>& arg)
  {
    std::string t;
    std::vector<T> v;
    is >> t;
    if (t != "[") GRAVIS_THROW3(gravis::Exception, "Unexpected token. A vector should start with [", t);
    while (is)
    {
      is >> t;
      if (t == "]") break;
      std::stringstream st(t);
      T tt;
      st >> tt;
      v.push_back(tt);
    }
    arg.resize(v.size());
    size_t h = arg.h;
    for (size_t j=0; j<h; ++j)
      arg[j] = v[j];
    return is;
  }
  /**
   * Write Variable size matrices to a stream
   **/
  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const tConstVectorView<T>& arg)
  {
    size_t h = arg.h;
    os << "[";
    for (size_t j=0; j<h; ++j)
      os << std::setw(8) << arg[j] << " ";
    os << " ]";
    return os;
  }
  /**
   * Write Variable size matrices to a stream
   **/
  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const tConstMatrixView<T>& arg)
  {
    size_t h = arg.h;
    size_t w = arg.w;
    if ((h>1) && (w>1))
    {
      os << "Matrix: " << h << "x" << w << std::endl;
      for (size_t i=0; i<h; ++i)
      {
        if (i==0)
        {
          os << "/";
        }
        else if (i==h-1)
        {
          os << "\\";
        }
        else
        {
          os << "|";
        }
        for (size_t j=0; j<w; ++j) os << std::setw(8) << arg(i,j) << " ";
        if (i==0)
        {
          os << "\\";
        }
        else if (i==h-1)
        {
          os << "/";
        }
        else
        {
          os << "|";
        }
        os << "\n";
      }
    }
    else if (w==1 && h>1)
    {
      os << "[";
      for (size_t j=0; j<h; ++j)
        os << std::setw(8) << arg[j] << " ";
      os << " ]^T";
    }
    else
    {
      os << "[";
      for (size_t j=0; j<w; ++j)
        os << std::setw(8) << arg[j] << " ";
      os << " ]";
    }
    return os;
  }

  /**
   * Write Variable size matrices to a stream
   **/
  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const tVarVector<T>& arg)
  {
    tConstVectorView<T> mv(arg);
    os << mv;
    return os;
  }
  /**
   * Write Variable size matrices to a stream
   **/
  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const tVarMatrix<T>& arg)
  {
    tConstMatrixView<T> mv(arg);
    os << mv;
    return os;
  }
}
#include "tVarMatrix_blas.h"

#endif
