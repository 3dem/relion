#ifndef __LIBGRAVIS_T_DEFAULT_VECTOR_H__
#define __LIBGRAVIS_T_DEFAULT_VECTOR_H__

#include <vector>
#include "tArray.h"

namespace gravis
{

  /**
   * Like a std::vector, but with a default value returned when accessing [-1].
   *
   * This situation is checked extremely efficiently by positioning the default
   * element at position [-1] in memory, so no check has to be done.
   *
   * This replacement for tVector does not offer reference counting. It makes
   * more sense to take a complete array structure and wrap it into a
   * boost::shared_ptr
   **/
  template <class T>
  class tDefaultVector
  {
    private:
      typedef typename std::vector<T> Vector;

      Vector data;
      T* data_ptr;


    public:
      typedef typename Vector::iterator                  iterator;
      typedef typename Vector::const_iterator            const_iterator;
      typedef typename Vector::reverse_iterator          reverse_iterator;
      typedef typename Vector::const_reverse_iterator    const_reverse_iterator;
      typedef typename Vector::reference                 reference;
      typedef typename Vector::const_reference           const_reference;

      /**
       * Create a new vector, optionally specifying a default value. If no default value is specified T() is used
       **/
      tDefaultVector(const size_t size=0, const T& def=T())
      {
        data.resize(size+1);
        data_ptr = &data[1];
        data[0] = def;
      }

      /**
       * Copy data from the other vector
       **/
      tDefaultVector(const tDefaultVector& other) : data(other.data), data_ptr(&data[1]) {};

      /**
       * Copy data from the other vector
       **/
      tDefaultVector(const tArray<T> &other) : data(other.size()+1), data_ptr(&data[1])
      {
        for (size_t i=0; i<other.size(); ++i)
          data[i+1] = other[i];
      };

      /**
       * Copy data from the other vector
       **/
      tDefaultVector(const std::vector<T> &other) : data(other.size()+1), data_ptr(&data[1])
      {
        for (size_t i=0; i<other.size(); ++i)
          data[i+1] = other[i];
      };

      /**
       * Exception save assignment operator
       **/
      void operator=(tDefaultVector& other)
      {
        Vector _data(other.data);
        std::swap(_data, data);
        data_ptr = &data[1];
      }

      inline size_t size() const
      {
        return data.size() - 1;
      }

      inline void push_back(const T& e)
      {
        data.push_back(e);
        data_ptr = &data[1];
      }
      inline iterator erase(iterator idx)
      {
        return data.erase(idx);
        data_ptr = &data[1];
      }
      inline iterator erase(iterator start, iterator end)
      {
        return data.erase(start, end);
        data_ptr = &data[1];
      }
      inline void erase(const size_t idx)
      {
        erase(begin()+idx);
        data_ptr = &data[1];
      }

      inline void resize(const size_t& size)
      {
        data.resize(size+1);
        data_ptr = &data[1];
      }

      /**
       * Without default checking
       **/
      inline reference operator[](const unsigned int& i)
      {
        return data_ptr[i];
      }
      /**
       * Without default checking
       **/
      inline const_reference operator[](const unsigned int& i) const
      {
        return data_ptr[i];
      }

      /**
       * With default checking
       **/
      inline reference operator[](const int& i)
      {
        return data_ptr[i];
      }
      /**
       * With default checking
       **/
      inline const_reference operator[](const int& i) const
      {
        return data_ptr[i];
      }

      /**
       * Without default checking
       **/
      inline reference operator[](const unsigned long& i)
      {
        return data_ptr[i];
      }
      /**
       * Without default checking
       **/
      inline const_reference operator[](const unsigned long& i) const
      {
        return data_ptr[i];
      }

      /**
       * With default checking
       **/
      inline reference operator[](const long& i)
      {
        return data_ptr[i];
      }
      /**
       * With default checking
       **/
      inline const_reference operator[](const long& i) const
      {
        return data_ptr[i];
      }

      inline iterator begin()
      {
        return(data.begin()+1);
      }
      inline iterator end()
      {
        return(data.end());
      }
      inline const_iterator begin() const
      {
        return(data.begin()++);
      }
      inline const_iterator end()   const
      {
        return(data.end());
      }

      inline void swap(tDefaultVector& other)
      {
        data.swap(other.data);
        std::swap(data_ptr, other.data_ptr);
      }

      inline void setDefault(const_reference def)
      {
        data[0] = def;
      };
      inline const_reference getDefault() const
      {
        return data[0];
      };
      inline reference getDefault()
      {
        return data[0];
      };

      inline void clear()
      {
        data.resize(1);
        data_ptr = &data[1];
      };
      inline void fill(const T& v)
      {
        for (size_t i=1; i<data.size(); ++i) data[i] = v;
      };

      /**
       * Reserve some additional space.
       **/
      inline void reserve(const size_t& space)
      {
        data.reserve(space+1);
        data_ptr = &data[1];
      }
  };


}

namespace std
{
  /// See tDefaultVector::swap().
  template<typename _Tp>
  inline void
  swap(gravis::tDefaultVector<_Tp>& __x, gravis::tDefaultVector<_Tp>& __y)
  {
    __x.swap(__y);
  }
}
#endif
