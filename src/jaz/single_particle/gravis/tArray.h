#ifndef __LIBGRAVIS_T_ARRAY_H__
#define __LIBGRAVIS_T_ARRAY_H__
/******************************************************************************
**        Title: tArray.h
**  Description: Implements a one dimensional array with reference counting.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <cassert>
#include <cstring>
#include <stddef.h>
#include <vector>
#include "private/tRefCPtr.h"


/*!
** \file tArray.h
*/

namespace gravis
{

  /*!
  ** \class tArray
  ** \brief Implements a one dimensional array with reference counting.
  */


  template <typename T>
  class tArray
  {
    public:
      typedef T value_type;

      tArray ();
      tArray (size_t);
      tArray (T* data, size_t nel, bool deleteData);
      tArray (const tArray&);
      tArray& operator=(const tArray&);
      tArray (const std::vector<T>&);
      ~tArray ();

      tArray<T> clone() const;
      tArray<T> safeClone() const;
      //! Deprecating this, as it is not standard. Use resize() instead.
      tArray<T>& setSize (size_t);
      //! Useful alias to setSize. At some point we should switch completely to the std library.
      tArray<T>& resize (size_t s)
      {
        return setSize(s);
      };

      void fill (T);
      void fill (T, size_t, size_t);

      size_t size () const;

      const T& operator[] (size_t) const;
      T& operator[] (size_t);

      const T* data () const;
      T* data ();

      bool operator==(const tArray& other) const;
      bool operator!=(const tArray& other) const;

      operator std::vector<T>() const
      {
        const tArray<T> &self = *this;
        const size_t l=size();
        std::vector<T> result(l);
        for (size_t i=0; i<l; ++i)
          result[i] = self[i];
        return result;
      }

    protected:
      void allocArray (size_t);

    protected:
      priv::tRefCPtr<T>  p_smp;
      T*                 p_data;
      size_t             length;
  };

  /*!
  ** \class tConstArray
  ** \brief Read-only wrapper for tArray.
  *
  * Since tArray is really a pointer, "const tArray" does protect the data,
  * but also protects the pointer! Assume I want a class that keeps a
  * pointer to data (= tArray), and needs only read access. We also want to
  * change the pointer once in a while.
  * \code
  * class X {
  *   // tArray<int> readOnly; // BAD! can manipulate data
  *   // const tArray<int> readOnly; // cannot manipulate data, but cannot change readOnly
  *   tConstArray<int> readOnly; // solution
  * public:
  *   void setArray(tConstArray<int> a) {
  *     readOnly = a;
  *   }
  * };
  * \endcode
  */
  template <typename T>
  class tConstArray
  {
    private:
      tArray<T> ta;
    public:
      tConstArray() {}
      tConstArray(tArray<T>& ta) : ta(ta) {}
      tArray<T> clone() const
      {
        return ta.clone();
      }
      tArray<T> safeClone() const
      {
        return ta.safeClone();
      }
      size_t size() const
      {
        return ta.size();
      }
      const T& operator[](size_t i) const
      {
        return ta[i];
      }
      const T* data() const
      {
        return ta.data();
      }
      bool operator==(const tArray<T>& other) const
      {
        return ta == other;
      }
      bool operator!=(const tArray<T>& other) const
      {
        return ta != other;
      }
      bool operator==(const tConstArray<T>& other) const
      {
        return ta == other.ta;
      }
      bool operator!=(const tConstArray<T>& other) const
      {
        return ta == other.ta;
      }
      const tConstArray& operator=(tArray<T>& ta)
      {
        this->ta = ta;
        return *this;
      }
      const tConstArray& operator=(const tConstArray<T>& other)
      {
        this->ta = other.ta;
        return *this;
      }
  };

  /*!
  ** \brief Default constructor
  */
  template <class T> inline tArray<T>::tArray () :
    p_smp (),
    p_data(),
    length()
  {
    this->allocArray(0);
  }


  /*!
  ** \brief Constructor.
  ** \param nel Number of elements to allocate for this tArray.
  */
  template <class T> inline tArray<T>::tArray (size_t nel) :
    p_smp (),
    p_data(),
    length()
  {
    this->allocArray(nel);
  }

  template <class T> inline tArray<T>::tArray (T* data, size_t nel, bool deleteData) :
    p_smp(),
    p_data(data),
    length(nel)
  {

    if (deleteData)
      p_smp  = priv::tRefCPtr<T>(p_data, priv::tRefCPtr<T>::ALLOC_ARRAY, 1);
    else
      p_smp  = priv::tRefCPtr<T>(p_data, priv::tRefCPtr<T>::ALLOC_ARRAY, 2);
  }


  /*!
  ** \brief Copy-constructor
  **
  ** The copy-constructor has reference-semantic, i.e. the managed data is not
  ** copied. Instead a new handle to the same data is created.
  **
  ** \param rhs The array to be copied
  */
  template <class T> inline tArray<T>::tArray (const tArray& rhs) :
    p_smp (rhs.p_smp),
    p_data(rhs.p_data),
    length(rhs.length)
  {
  }

  /*!
  ** \brief Assignment
  **
  ** The assignment has reference-semantic, i.e. the managed data is not
  ** copied. Instead a new handle to the same data is created.
  **
  ** \param rhs  The array to be assigned
  */
  template <class T> inline tArray<T> &tArray<T>::operator=(const tArray& rhs)
  {
    p_smp  = rhs.p_smp;
    p_data = rhs.p_data;
    length = rhs.length;
    return *this;
  }

  /*!
  ** \brief Construct from std vector
  **
  ** \param rhs The std vector from which the data is copied. This construction does not create a reference, but actually copies the data.
  */
  template <class T> inline tArray<T>::tArray (const std::vector<T>& rhs)
  {
    this->allocArray(rhs.size());
    for (size_t i=0; i<rhs.size(); ++i)
      (*this)[i] = rhs[i];
  }


  /*!
  ** \brief Destructor.
  **
  ** Destroy the object. The managed data is *only* deleted if no other
  ** tArray object holds a reference to it.
  */
  template <class T> inline tArray<T>::~tArray ()
  {
  }

  /*!
  ** \brief Create a deep-copy of managed data.

  ** \return A new tArray<T> object.
  **
  ** Use this version of clone unless your datatype is simple
  ** (e.g. tVector, size_t, Tuple2...)
  */
  template <class T> inline tArray<T> tArray<T>::safeClone() const
  {
    tArray lhs(length);
    for (size_t i=0; i<length; i++)
      lhs.p_data[i] = p_data[i];
    return lhs;
  }

  /*!
  ** \brief Create a deep-copy of managed data.
  **
  ** \return A new tArray<T> object.
  **
  ** \warning This method creates a byte-wise copy of the managed data. When
  ** applied to compound types (e.g. T=std::vector<int>) or reference counted
  ** types like std::string it will create crashes use save_clone() unless your
  datatype is simple.
  */
  template <class T> inline tArray<T> tArray<T>::clone () const
  {
    tArray lhs(length);
    memcpy(lhs.p_data, p_data, length*sizeof(T));
    return lhs;
  }


  /*!
  ** \brief Fill array with constant value.
  **
  ** \param value Value to fill with.
  */
  template <class T> inline void tArray<T>::fill (T value)
  {
    const T* end_ptr = p_data + length;
    for (T* t_ptr=p_data; t_ptr<end_ptr; ++t_ptr) *t_ptr = value;
  }

  /*!
  ** \brief Fill array with constant value.
  **
  ** \param value Value to fill with.
  ** \param from  Starting position (defaults to first element).
  ** \param to    Last position (defaults to last element).
  */
  template <class T> inline void tArray<T>::fill (T value, size_t from, size_t to)
  {
    if (from >= length) from = length;
    if (to   >= length) to   = length;

    T* end_ptr = p_data + to;
    for (T* t_ptr=p_data+from; t_ptr<end_ptr; ++t_ptr) *t_ptr = value;
  }

  /*!
  ** \brief Resize array.
  ** \todo rename to resize() (the name in all containers, including valarray, which has even the same semantics)
  ** This method changes the number of T-elements managed by the object.
  **
  ** \param nel Number of T elements in resized array.
  ** \return Reference to resized array object.
  ** \warning The data managed by the array object is not copied.
  */
  template <class T> inline tArray<T>& tArray<T>::setSize (size_t nel)
  {
    this->allocArray(nel);
    return *this;
  }


  /*!
  ** \brief Get number of elements.
  **
  ** \return Number of T-elements in array.
  */
  template <class T> inline size_t tArray<T>::size () const
  {
    return length;
  }


  /*!
  ** \brief Access i-th element.
  ** \param i Index into array.
  ** \return const-Reference to i-th element.
  */
  template <class T> inline const T& tArray<T>::operator[] (size_t i) const
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert( i < length );
#endif
    return p_data[i];
  }


  /*!
  ** \brief Access i-th element.
  ** \param i Index into array.
  ** \return Reference to i-th element.
  */
  template <class T> inline T& tArray<T>::operator[] (size_t i)
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert( i < length );
#endif
    return p_data[i];
  }

  /*!
  ** \brief Perform element-by-element comparison.
  */
  template <class T> inline bool tArray<T>::operator==(const tArray<T>& other) const
  {
    return !(*this != other);
  }

  /*!
  ** \brief Perform element-by-element comparison.
  */
  template <class T> inline bool tArray<T>::operator!=(const tArray<T>& other) const
  {
    if (p_data == other.p_data) return false;
    else if (length != other.length) return true;
    else
    {
      for (size_t i = 0; i < length; i++)
      {
        if (p_data[i] != other.p_data[i]) return true;
      }
    }
    return false;
  }

  /*!
  ** \brief Get pointer to managed data.
  ** \return const-Pointer to first element of managed data.
  */
  template <class T> inline const T* tArray<T>::data () const
  {
    return p_data;
  }


  /*!
  ** \brief Get pointer to managed data.
  ** \return Pointer to first element of managed data.
  */
  template <class T> inline T* tArray<T>::data ()
  {
    return p_data;
  }


  template <class T> inline void tArray<T>::allocArray (size_t nel)
  {
    if (nel <= 0)
    {
      p_data = 0;
      length = 0;
      p_smp  = priv::tRefCPtr<T>(p_data, priv::tRefCPtr<T>::ALLOC_ARRAY);
    }
    else
    {
      // ATTENTION! Bug: this leaks!! ... delete old memory
      // Sandro Schoenborn, 2013-04-09, sandro.schoenborn@unibas.ch
      // Tobias Maier, 2013-04-09, tobias.maier@unibas.ch
      p_data = new T[nel];
      length = nel;
      p_smp  = priv::tRefCPtr<T>(p_data, priv::tRefCPtr<T>::ALLOC_ARRAY);
    }
  }

} /* Close namespace "gravis" */

#endif
