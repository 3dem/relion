#ifndef __LIBGRAVIS_NTUPLE_H__
#define __LIBGRAVIS_NTUPLE_H__
/******************************************************************************
 **        Title: NNTuple.h
 **  Description: N-NTuples (templated)
 **
 **       Author: Brian Amberg 2006
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

namespace gravis
{

  /*! \brief N-NTuple, typically used with Ieger datatypes for multi-index. */
  template <class I, size_t N>
  class NTuple
  {
      I m[N];

    public:
      /*! \brief Construct a NTuple with entries of -1. */
      NTuple()
      {
        for (size_t i=0; i<N; ++i) (*this)[i] = -1;
      }
      const I& operator[](size_t i) const
      {
        return m[i];
      }
      I& operator[](size_t i)
      {
        return m[i];
      }

      /*! \brief Offset all non-negative entries. */
      void offset(I o)
      {
        for (size_t i=0; i<N; ++i)
          if ((*this)[i] >= I(0))
            (*this)[i] += o;
      }

      /*! \brief Whether all entries are non-negative. */
      bool allValid() const
      {
        bool r = (*this)[0] >= I(0);
        for (size_t i=1; i<N; ++i)
          r = r && ((*this)[i] >= I(0));
        return r;
      }

      //! Lexical Ordering for NTuples
      inline bool operator==(const NTuple& o) const
      {
        bool r = (*this)[0] == o[0];
        for (size_t i=1; i<N; ++i)
          r = r && ((*this)[i] == o[i]);
        return r;
      }

      //! Lexical Ordering for NTuples
      inline bool operator!=(const NTuple& o) const
      {
        return !(*this == o);
      }

      //! Lexical Ordering for NTuples
      inline bool operator<(const NTuple& o) const
      {
        for (size_t i=0; i<N; ++i)
          if ((*this)[i] < o[i])
            return true;
        return false;
      }

      //! Lexical Ordering for NTuples
      inline bool operator>(const NTuple& o) const
      {
        return (*this != o) && !(*this < o);
      }

      //! Lexical Ordering for NTuples
      inline bool operator<=(const NTuple& o) const
      {
        return (*this < o) || (*this == o);
      }

      //! Lexical Ordering for NTuples
      inline bool operator>=(const NTuple& o) const
      {
        return (*this > o) || (*this == o);
      }
  }; // class NTuple


  template <class I, unsigned int N>
  inline
  std::ostream& operator<< (std::ostream& os, const NTuple<I, N>& arg)
  {
    os << "[";
    for (int i=0; i<N-1; ++i) os << arg[i] << ", ";
    os << arg[N-1] << "]";
    return os;
  }

  template <class I> NTuple<I, 1> nTuple(const I& c0)
  {
    NTuple<I, 1> r;
    r[0] = c0;
    return r;
  };
  template <class I> NTuple<I, 2> nTuple(const I& c0, const I& c1)
  {
    NTuple<I, 2> r;
    r[0] = c0;
    r[1] = c1;
    return r;
  };
  template <class I> NTuple<I, 3> nTuple(const I& c0, const I& c1, const I& c2)
  {
    NTuple<I, 3> r;
    r[0] = c0;
    r[1] = c1;
    r[2] = c2;
    return r;
  };
  template <class I> NTuple<I, 3> nTuple(const I& c0, const I& c1, const I& c2, const I& c3)
  {
    NTuple<I, 3> r;
    r[0] = c0;
    r[1] = c1;
    r[2] = c2;
    r[3] = c3;
    return r;
  };

  typedef NTuple<int, 1> I1Tuple;
  typedef NTuple<int, 2> I2Tuple;
  typedef NTuple<int, 3> I3Tuple;
  typedef NTuple<int, 4> I4Tuple;
} // namespace gravis

#endif
