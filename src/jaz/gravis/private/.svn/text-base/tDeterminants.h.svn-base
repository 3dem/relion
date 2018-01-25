/******************************************************************************
 **        Title: tDeterminants.h
 **  Description: Templated functions for 2D, 3D and 4D determinants.
 **
 **       Author: Michael Keller, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#ifndef _TDETERMINANTS_H_
#define _TDETERMINANTS_H_

namespace gravis
{

  template <class T> inline
  T det2x2(T a1, T a2,
           T b1, T b2)
  {
    return
      + a1 * b2
      - b1 * a2;
  }

  template <class T> inline
  T det3x3(T a1, T a2, T a3,
           T b1, T b2, T b3,
           T c1, T c2, T c3)
  {
    return
      + a1 * det2x2(b2, b3, c2, c3)
      - b1 * det2x2(a2, a3, c2, c3)
      + c1 * det2x2(a2, a3, b2, b3);
  }

  template <class T> inline
  T det4x4(T a1, T a2, T a3, T a4,
           T b1, T b2, T b3, T b4,
           T c1, T c2, T c3, T c4,
           T d1, T d2, T d3, T d4)
  {
    return
      + a1 * det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
      - b1 * det3x3(c2, c3, c4, d2, d3, d4, a2, a3, a4)
      + c1 * det3x3(d2, d3, d4, a2, a3, a4, b2, b3, b4)
      - d1 * det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4)
      - a2 * det3x3(b3, b4, b1, c3, c4, c1, d3, d4, d1)
      + b2 * det3x3(c3, c4, c1, d3, d4, d1, a3, a4, a1)
      - c2 * det3x3(d3, d4, d1, a3, a4, a1, b3, b4, b1)
      + d2 * det3x3(a3, a4, a1, b3, b4, b1, c3, c4, c1)
      + a3 * det3x3(b4, b1, b2, c4, c1, c2, d4, d1, d2)
      - b3 * det3x3(c4, c1, c2, d4, d1, d2, a4, a1, a2)
      + c3 * det3x3(d4, d1, d2, a4, a1, a2, b4, b1, b2)
      - d3 * det3x3(a4, a1, a2, b4, b1, b2, c4, c1, c2)
      - a4 * det3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3)
      + b4 * det3x3(c1, c2, c3, d1, d2, d3, a1, a2, a3)
      - c4 * det3x3(d1, d2, d3, a1, a2, a3, b1, b2, b3)
      + d4 * det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
  }

}

#endif
