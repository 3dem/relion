#ifndef __LIBGRAVIS_TUPLE_H__
#define __LIBGRAVIS_TUPLE_H__
/******************************************************************************
 **        Title: Tuple.h
 **  Description: Tuple2 and Tuple3 (tuples of int)
 **
 **       Author: Michael Keller, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace gravis
{

  /*! \brief Tuple of 2 integers, typically used for multi-index. */
  class Tuple2
  {

    public:
      int c0, c1;
      /*! \brief Construct a Tuple2 with entries of -1. */
      Tuple2() : c0(-1), c1(-1) {}
      Tuple2(int a, int b) : c0(a), c1(b) {}
      int operator[] (int i) const
      {
        return *(&c0 + i);
      }
      int& operator[] (int i)
      {
        return *(&c0 + i);
      }
      unsigned int size() const
      {
        return 2;
      }

      /*! \brief Offset all non-negative entries. */
      void offset(int o)
      {
        if (c0 >= 0) c0 += o;
        if (c1 >= 0) c1 += o;
      }
      /*! \brief Whether all entries are non-negative. */
      bool allValid() const
      {
        return c0 >= 0 && c1 >= 0;
      }

      //! Lexical Ordering for Tuples
      inline bool operator==(const Tuple2& o) const
      {
        return ((c0 == o.c0) && (c1 == o.c1));
      }

      //! Lexical Ordering for Tuples
      inline bool operator!=(const Tuple2& o) const
      {
        return !(*this == o);
      }

      //! Lexical Ordering for Tuples
      inline bool operator<(const Tuple2& o) const
      {
        return ((c0 < o.c0) || ((c0 == o.c0) && (c1 < o.c1)));
      }

      //! Lexical Ordering for Tuples
      inline bool operator>(const Tuple2& o) const
      {
        return (*this != o) && !(*this < o);
      }

      //! Lexical Ordering for Tuples
      inline bool operator<=(const Tuple2& o) const
      {
        return (*this < o) || (*this == o);
      }

      //! Lexical Ordering for Tuples
      inline bool operator>=(const Tuple2& o) const
      {
        return (*this > o) || (*this == o);
      }
  }; // class Tuple2

  /*! \brief Tuple of three integers, typically used for multi-index. */
  class Tuple3
  {

    public:
      int c0, c1, c2;
      /*! \brief Construct a Tuple3 with entries of -1. */
      Tuple3() : c0(-1), c1(-1), c2(-1) {}
      Tuple3(int a, int b, int c) : c0(a), c1(b), c2(c) {}
      int operator[] (int i) const
      {
        return *(&c0 + i);
      }
      int& operator[] (int i)
      {
        return *(&c0 + i);
      }
      unsigned int size() const
      {
        return 3;
      }

      /*! brief Offset all non-negative entries. */
      void offset(int o)
      {
        if (c0 >= 0) c0 += o;
        if (c1 >= 0) c1 += o;
        if (c2 >= 0) c2 += o;
      }
      /*! \brief Whether all entries are non-negative. */
      bool allValid() const
      {
        return c0 >= 0 && c1 >= 0 && c2 >= 0;
      }

      //! Lexical Ordering for Tuples
      inline bool operator==(const Tuple3& o) const
      {
        return ((c0 == o.c0) && (c1 == o.c1) && (c2 == o.c2));
      }

      //! Lexical Ordering for Tuples
      inline bool operator!=(const Tuple3& o) const
      {
        return !(*this == o);
      }

      //! Lexical Ordering for Tuples
      inline bool operator<(const Tuple3& o) const
      {
        return ((c0 < o.c0) || ((c0 == o.c0) && (c1 < o.c1)) || ((c0 == o.c0) && (c1 == o.c1) && (c2 < o.c2)));
      }

      //! Lexical Ordering for Tuples
      inline bool operator>(const Tuple3& o) const
      {
        return (*this != o) && !(*this < o);
      }

      //! Lexical Ordering for Tuples
      inline bool operator<=(const Tuple3& o) const
      {
        return (*this < o) || (*this == o);
      }

      //! Lexical Ordering for Tuples
      inline bool operator>=(const Tuple3& o) const
      {
        return (*this > o) || (*this == o);
      }
  }; // class Tuple3

  inline
  std::ostream& operator<< (std::ostream& os, const Tuple3& arg)
  {
    os << "[" << arg.c0 << ", " << arg.c1 << ", " << arg.c2 << "]";
    return os;
  }

  inline
  std::ostream& operator<< (std::ostream& os, const Tuple2& arg)
  {
    os << "[" << arg.c0 << ", " << arg.c1 << "]";
    return os;
  }
  
  // Inverse of operator<<
  inline
  std::istream& operator>> (std::istream& is, Tuple3& arg)
  {
    char c = ' ';
    is >> c;
    if (is.eof())
      return is;
    if (c != '[')
      throw std::runtime_error("Tuple should start with an opening [");
    std::stringstream values;
    int v = 0;
    while ((is >> c) && (c != ']'))
    {
      if (c == ',')
      {
        v++;
        if (v >= 3)
          throw std::runtime_error("Tuple3 contains more than three elements");
        values << " ";
      }
      else if (c != ' ')
        values << c;
    }
    if (c != ']')
    {
      throw std::runtime_error("Tuple3 should end with a ]");
    }
    values >> arg.c0 >> arg.c1 >> arg.c2;
    return is;
  }

    // Inverse of operator<<
  inline
  std::istream& operator>> (std::istream& is, Tuple2& arg)
  {
    char c = ' ';
    is >> c;
    if (is.eof())
      return is;
    if (c != '[')
      throw std::runtime_error("Tuple should start with an opening [");
    std::stringstream values;
    int v = 0;
    while ((is >> c) && (c != ']'))
    {
      if (c == ',')
      {
        v++;
        if (v >= 2)
          throw std::runtime_error("Tuple2 contains more than three elements");
        values << " ";
      }
      else if (c != ' ')
        values << c;
    }
    if (c != ']')
    {
      throw std::runtime_error("Tuple2 should end with a ]");
    }
    values >> arg.c0 >> arg.c1;
    return is;
  }

} // namespace gravis

#endif
