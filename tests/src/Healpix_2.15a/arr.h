/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file src/Healpix_2.15a/cxxsupport/arr.h
 *  Various high-performance array classes used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ARR_H
#define PLANCK_ARR_H

#include "src/Healpix_2.15a/cxxutils.h"
#include <algorithm>

/*! \defgroup arraygroup Array classes */
/*! \{ */

/*! An array whose size is known at compile time. Very useful for storing
    small arrays on the stack, without need for \a new and \a delete(). */
template <typename T, unsigned int sz> class fix_arr
  {
  private:
    T d[sz];

  public:
    /*! Returns the size of the array. */
    long size() const { return sz; }

    /*! Returns a reference to element \a #n */
    template<typename T2> T &operator[] (T2 n) {return d[n];}
    /*! Returns a constant reference to element \a #n */
    template<typename T2> const T &operator[] (T2 n) const {return d[n];}
  };


/*! One-dimensional array type. */
template <typename T> class arr
  {
  private:
    long s;
    T *d;
    bool own;

#if defined(PLANCK_CHECKS)
    void check_range(long n) const
      {
      if ((n<0) || (n>=s)) throw Message_error
        ("arr: index "+dataToString(n)+" is out of range. Max index is "
         +dataToString(s-1));
      }
#endif

    void reset()
      { s=0; d=0; own=true; }

  public:
    /*! Creates a zero-sized array. */
    arr() : s(0), d(0), own(true) {}
    /*! Creates an array with \a sz entries. */
    arr(long sz) : s(sz), d (s>0 ? new T[s] : 0), own(true) {}
    /*! Creates an array with \a sz entries, and initializes them with
        \a inival. */
    arr(long sz, const T &inival) : s(sz), d (s>0 ? new T[s] : 0), own(true)
      { fill(inival); }
    /*! Creates an array with \a sz entries, which uses the memory pointed
        to by \a ptr.
        \note \a ptr will <i>not</i> be deallocated by the destructor.
        \warning Only use this if you REALLY know what you are doing.
        In particular, this is only safely usable if
          <ul>
          <li>\a T is a POD type</li>
          <li>\a ptr survives during the lifetime of the array object</li>
          <li>\a ptr is not subject to garbage collection</li>
          </ul>
        Other restrictions may apply. You have been warned. */
    arr (T *ptr, long sz): s(sz), d(ptr), own(false) {}
    /*! Creates an array which is a copy of \a orig. The data in \a orig
        is duplicated. */
    arr (const arr &orig): s(orig.s), d (s>0 ? new T[s] : 0), own(true)
      { for (long m=0; m<s; ++m) d[m] = orig.d[m]; }
    /*! Frees the memory allocated by the object. */
    ~arr() { if (own) delete[] d; }

    /*! Returns the current array size. */
    long size() const { return s; }

    /*! Allocates space for \a sz elements. The content of the array is
        undefined on exit. \a sz can be 0. If \a sz is the
        same as the current size, no reallocation is performed. */
    void alloc (long sz)
      {
      if (sz==s) return;
      if (own) delete[] d;
      s = sz;
      d = s>0 ? new T[sz] : 0;
      own = true;
      }
    /*! Deallocates the memory held by the array, and sets the array size
        to 0. */
    void dealloc () {if (own) delete[] d; reset();}

    /*! Writes \a val into every element of the array. */
    void fill (const T &val)
      { for (long m=0; m<s; ++m) d[m]=val; }

    /*! Changes the array to be a copy of \a orig. */
    arr &operator= (const arr &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s);
      for (long m=0; m<s; ++m) d[m] = orig.d[m];
      return *this;
      }

#if defined (PLANCK_CHECKS)
    template<typename T2> T &operator[] (T2 n) {check_range(n); return d[n];}
    template<typename T2> const T &operator[] (T2 n) const
      {check_range(n); return d[n];}
#else
    /*! Returns a reference to element \a #n */
    template<typename T2> T &operator[] (T2 n) {return d[n];}
    /*! Returns a constant reference to element \a #n */
    template<typename T2> const T &operator[] (T2 n) const {return d[n];}
#endif

    T *begin() { return d; }
    T *end() { return d+s; }

    /*! Sorts the elements in the array, in ascending order. */
    void sort()
      { std::sort (d,d+s); }

    /*! Returns the minimum and maximum entry in \a minv and \a maxv,
        respectively. Throws an exception if the array is zero-sized. */
    void minmax (T &minv, T &maxv) const
      {
      planck_assert(s>0,"trying to find min and max of a zero-sized array");
      minv=maxv=d[0];
      for (int m=1; m<s; ++m)
        {
        if (d[m]<minv) minv=d[m];
        else if (d[m]>maxv) maxv=d[m];
        }
      }

    /*! Assigns the contents and size of \a other to the array. On exit,
        \a other is yero-sized. */
    void transfer (arr &other)
      { if (own) delete[] d; d=other.d; s=other.s; own=other.own; other.reset(); }
    /*! Swaps contents and size with \a other. */
    void swap (arr &other)
      { std::swap(d,other.d); std::swap(s,other.s); std::swap(own,other.own);}
  };

/*! Two-dimensional array type. The storage ordering is the same as in C.
    An entry is located by address arithmetic, not by double dereferencing.
    The indices start at zero. */
template <typename T> class arr2
  {
  private:
    long s1, s2;
    arr<T> d;

#if defined (PLANCK_CHECKS)
    void check_range(long n) const
      {
      if ((n<0) || (n>=s1)) throw Message_error
        ("arr2: index "+dataToString(n)+" is out of range. Max index is "
         +dataToString(s1-1));
      }
#endif

  public:
    /*! Creates a zero-sized array. */
    arr2() : s1(0), s2(0) {}
    /*! Creates an array with the dimensions \a sz1 and \a sz2. */
    arr2(long sz1, long sz2)
      : s1(sz1), s2(sz2), d(s1*s2) {}
    /*! Creates the array as a copy of \a orig. */
    arr2(const arr2 &orig)
      : s1(orig.s1), s2(orig.s2), d(orig.d) {}
    /*! Frees the memory associated with the array. */
    ~arr2() {}

    /*! Returns the first array dimension. */
    long size1() const { return s1; }
    /*! Returns the second array dimension. */
    long size2() const { return s2; }
    /*! Returns the total array size, i.e. the product of both dimensions. */
    long size () const { return s1*s2; }

    /*! Allocates space for an array with \a sz1*sz2 elements.
        The content of the array is undefined on exit.
        \a sz1 or \a sz2 can be 0. If \a sz1*sz2 is the same as the
        currently allocated space, no reallocation is performed. */
    void alloc (long sz1, long sz2)
      {
      if (sz1*sz2 != d.size())
        d.alloc(sz1*sz2);
      s1=sz1; s2=sz2;
      }
    /*! Allocates space for an array with \a sz1*sz2 elements.
        The content of the array is undefined on exit.
        \a sz1 or \a sz2 can be 0. If \a sz1*sz2 is smaller than the
        currently allocated space, no reallocation is performed. */
    void fast_alloc (long sz1, long sz2)
      {
      if (sz1*sz2<=d.size())
        { s1=sz1; s2=sz2; }
      else
        alloc(sz1,sz2);
      }
    /*! Deallocates the space and makes the array zero-sized. */
    void dealloc () {d.dealloc(); s1=0; s2=0;}

    /*! Sets all array elements to \a val. */
    void fill (const T &val)
      { d.fill(val); }

    /*! Changes the array to be a copy of \a orig. */
    arr2 &operator= (const arr2 &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s1, orig.s2);
      d = orig.d;
      return *this;
      }

#if defined (PLANCK_CHECKS)
    template<typename T2> T *operator[] (T2 n)
      {check_range(n);return &d[n*s2];}
    template<typename T2> const T *operator[] (T2 n) const
      {check_range(n);return &d[n*s2];}
#else
    /*! Returns a pointer to the beginning of slice \a #n. */
    template<typename T2> T *operator[] (T2 n) {return &d[n*s2];}
    /*! Returns a constant pointer to the beginning of slice \a #n. */
    template<typename T2> const T *operator[] (T2 n) const {return &d[n*s2];}
#endif

    /*! Returns the minimum and maximum entry in \a minv and \a maxv,
        respectively. Throws an exception if the array is zero-sized. */
    void minmax (T &minv, T &maxv) const
      {
      planck_assert(s1*s2>0,
        "trying to find min and max of a zero-sized array");
      minv=maxv=d[0];
      for (int m=1; m<s1*s2; ++m)
        {
        if (d[m]<minv) minv=d[m];
        if (d[m]>maxv) maxv=d[m];
        }
      }

    /*! Swaps contents and sizes with \a other. */
    void swap (arr2 &other)
      {
      d.swap(other.d);
      std::swap(s1,other.s1);
      std::swap(s2,other.s2);
      }
  };

/*! Two-dimensional array type. An entry is located by double dereferencing,
    i.e. via an array of pointers. The indices start at zero. */
template <typename T> class arr2b
  {
  private:
    long s1, s2;
    arr<T> d;
    arr<T *> d1;

#if defined (PLANCK_CHECKS)
    void check_range(long n) const
      {
      if ((n<0) || (n>=s1)) throw Message_error
        ("arr: index "+dataToString(n)+" is out of range. Max index is "
         +dataToString(s1-1));
      }
#endif

    void fill_d1()
      { for (long m=0; m<s1; ++m) d1[m] = &d[m*s2]; }

  public:
    /*! Creates a zero-sized array. */
    arr2b() : s1(0), s2(0), d(0), d1(0) {}
    /*! Creates an array with the dimensions \a sz1 and \a sz2. */
    arr2b(long sz1, long sz2)
      : s1(sz1), s2(sz2), d(s1*s2), d1(s1)
      { fill_d1(); }
    /*! Creates the array as a copy of \a orig. */
    arr2b(const arr2b &orig)
      : s1(orig.s1), s2(orig.s2), d(orig.d), d1(s1)
      { fill_d1(); }
    /*! Frees the memory associated with the array. */
    ~arr2b() {}

    /*! Returns the first array dimension. */
    long size1() const { return s1; }
    /*! Returns the second array dimension. */
    long size2() const { return s2; }
    /*! Returns the total array size, i.e. the product of both dimensions. */
    long size () const { return s1*s2; }

    /*! Allocates space for an array with \a sz1*sz2 elements.
        The content of the array is undefined on exit. */
    void alloc (long sz1, long sz2)
      {
      if ((s1==sz1) && (s2==sz2)) return;
      s1=sz1; s2=sz2;
      d.alloc(s1*s2);
      d1.alloc(s1);
      fill_d1();
      }
    /*! Deallocates the space and makes the array zero-sized. */
    void dealloc () {d.dealloc(); d1.dealloc(); s1=0; s2=0;}

    /*! Sets all array elements to \a val. */
    void fill (const T &val)
      { d.fill(val); }

    /*! Changes the array to be a copy of \a orig. */
    arr2b &operator= (const arr2b &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s1, orig.s2);
      for (long m=0; m<s1*s2; ++m) d[m] = orig.d[m];
      return *this;
      }

#if defined (PLANCK_CHECKS)
    template<typename T2> T *operator[] (T2 n) {check_range(n); return d1[n];}
    template<typename T2> const T *operator[] (T2 n) const
      {check_range(n); return d1[n];}
#else
    /*! Returns a pointer to the beginning of slice \a #n. */
    template<typename T2> T *operator[] (T2 n) {return d1[n];}
    /*! Returns a constant pointer to the beginning of slice \a #n. */
    template<typename T2> const T *operator[] (T2 n) const {return d1[n];}
#endif
    /*! Returns a pointer to the beginning of the pointer array. */
    T **p0() {return &d1[0];}
  };


/*! Three-dimensional array type. The storage ordering is the same as in C.
    An entry is located by address arithmetic, not by multiple dereferencing.
    The indices start at zero. */
template <typename T> class arr3
  {
  private:
    long s1, s2, s3, s2s3;
    arr<T> d;

  public:
    /*! Creates a zero-sized array. */
    arr3() : s1(0), s2(0), s3(0), s2s3(0), d(0) {}
    /*! Creates an array with the dimensions \a sz1, \a sz2 and \a sz3. */
    arr3(long sz1, long sz2, long sz3)
      : s1(sz1), s2(sz2), s3(sz3), s2s3(s2*s3), d(s1*s2*s3) {}
    /*! Creates the array as a copy of \a orig. */
    arr3(const arr3 &orig)
      : s1(orig.s1), s2(orig.s2), s3(orig.s3), s2s3(orig.s2s3), d(orig.d) {}
    /*! Frees the memory associated with the array. */
    ~arr3() {}

    /*! Returns the first array dimension. */
    long size1() const { return s1; }
    /*! Returns the second array dimension. */
    long size2() const { return s2; }
    /*! Returns the third array dimension. */
    long size3() const { return s3; }
    /*! Returns the total array size, i.e. the product of all dimensions. */
    long size () const { return s1*s2*s3; }

    /*! Allocates space for an array with \a sz1*sz2*sz3 elements.
        The content of the array is undefined on exit. */
    void alloc (long sz1, long sz2, long sz3)
      {
      d.alloc(sz1*sz2*sz3);
      s1=sz1; s2=sz2; s3=sz3; s2s3=s2*s3;
      }
    /*! Deallocates the space and makes the array zero-sized. */
    void dealloc () {d.dealloc(); s1=0; s2=0; s3=0; s2s3=0;}

    /*! Sets all array elements to \a val. */
    void fill (const T &val)
      { d.fill(val); }

    /*! Changes the array to be a copy of \a orig. */
    arr3 &operator= (const arr3 &orig)
      {
      if (this==&orig) return *this;
      alloc (orig.s1, orig.s2, orig.s3);
      d = orig.d;
      return *this;
      }

    /*! Returns a reference to the element with the indices
        \a n1, \a n2 and \a n3. */
    template<typename T2> T &operator() (T2 n1, T2 n2, T2 n3)
      {return d[n1*s2s3 + n2*s3 + n3];}
    /*! Returns a constant reference to the element with the indices
        \a n1, \a n2 and \a n3. */
    template<typename T2> const T &operator() (T2 n1, T2 n2, T2 n3) const
      {return d[n1*s2s3 + n2*s3 + n3];}

    /*! Swaps contents and sizes with \a other. */
    void swap (arr3 &other)
      {
      d.swap(other.d);
      std::swap(s1,other.s1);
      std::swap(s2,other.s2);
      std::swap(s3,other.s3);
      std::swap(s2s3,other.s2s3);
      }
  };

/*! \} */

#endif
