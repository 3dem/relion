#ifndef __LIBGRAVIS_T_QUATERNION_H__
#define __LIBGRAVIS_T_QUATERNION_H__
/******************************************************************************
 **        Title: tQuaternion.h
 **  Description: Represents a quaternion useful for rotation and scaling
 **
 **       Author: Reinhard Knothe
 **               Brian Amberg
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#include "t4Matrix.h"

namespace gravis
{

  /**
   * A tQuaternion class useful for rotation+scaling operations
   **/
  template<class T>
  class tQuaternion
  {

    public:
      T s,v1,v2,v3;

      tQuaternion(const T& s=0, const T& v1=1, const T& v2=0, const T& v3=0) : s(s), v1(v1), v2(v2), v3(v3) {}
      tQuaternion(const tQuaternion& q) : s(q.s), v1(q.v1), v2(q.v2), v3(q.v3) {}
      tQuaternion(const T* q) : s(q[0]), v1(q[1]), v2(q[2]), v3(q[3]) {}

      tQuaternion(const T& phi, const gravis::t3Vector<T> &axis) :
        s(cos(phi/T(2))), v1(axis.x* sin(phi/T(2))), v2(axis.y* sin(phi/T(2))), v3(axis.z* sin(phi/T(2))) {}

      bool operator==(const tQuaternion& q) const
      {
        return (s==q.s) && (v1==q.v1) && (v2==q.v2) && (v3==q.v3);
      }
	  
      bool operator!=(const tQuaternion& q) const
      {
        return !(*this == q);
      }

      tQuaternion operator + (const tQuaternion& q) const
      {
        return tQuaternion(q.s+s, q.v1+v1, q.v2+v2, q.v3+v3 );
      }
	  
      tQuaternion operator - (const tQuaternion& q) const
      {
        return tQuaternion(s-q.s, v1-q.v1, v2-q.v2, v3-q.v3 );
      }
	  
      tQuaternion operator * (const tQuaternion& q) const
      {
        t3Vector<T> v00 (v1,v2,v3);
        t3Vector<T> v10 (q.v1,q.v2,q.v3);

        T s2 = s * q.s - dot(v00,v10);
        t3Vector<T> v2 = cross(v00,v10);
        t3Vector<T> v3 = v10;
        v3 *= s;
        t3Vector<T> v4 = v00;
        v4 *= q.s;
        t3Vector<T> v5 = v2+v3+v4;

        return tQuaternion(s2, v5.x, v5.y, v5.z);
      }

      /**
       * Norm
       **/
      T length() const
      {
        return sqrt(s*s + v1*v1 + v2*v2 + v3*v3);
      }

      /**
       * Inplace normalization
       **/
      void normalize()
      {
        T l = length();
        s /= l;
        v1 /= l;
        v2 /= l;
        v3 /= l;
      }

      t3Matrix<T> getMatrix3() const
      {
        return t3Matrix<T>(
                 T(1)-T(2)*(v2*v2 + v3*v3),  T(2)*(v1*v2 - v3*s),            T(2)*(v3*v1 + v2*s),
                 T(2) * (v1*v2 + v3*s),      T(1) - T(2) * (v3*v3 + v1*v1),  T(2) * (v2*v3 - v1*s),
                 T(2) * (v3*v1 - v2*s),      T(2) * (v2*v3 + v1*s),          T(1) - T(2) * (v2*v2 + v1*v1));
      }

      t4Matrix<T> getMatrix4() const
      {
        return t4Matrix<T>(
                 T(1)-T(2)*(v2*v2 + v3*v3),  T(2)*(v1*v2 - v3*s),            T(2)*(v3*v1 + v2*s),            T(0),
                 T(2) * (v1*v2 + v3*s),      T(1) - T(2) * (v3*v3 + v1*v1),  T(2) * (v2*v3 - v1*s),          T(0),
                 T(2) * (v3*v1 - v2*s),      T(2) * (v2*v3 + v1*s),          T(1) - T(2) * (v2*v2 + v1*v1),  T(0),
                 T(0),                       T(0),                           T(0),                           T(1)  );
      }
  };

  template <class T>
  inline
  std::ostream& operator<< (std::ostream& os, const tQuaternion<T>& arg)
  {
    os << "[" << arg.s << "; " << arg.v1 << ", " << arg.v2 << ", " << arg.v3 << "]";
    return os;
  }

  typedef gravis::tQuaternion<float>  fQuaternion;
  typedef gravis::tQuaternion<double> dQuaternion;
}
#endif
