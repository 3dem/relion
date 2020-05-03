/*************************************************************************//*!
*        Title: tImage/operators.h
*  Description: Implements operators on tImages
*
*       Author: Brian Amberg, 2006-2007
*               Computer Science Department, University Basel (CH)
****************************************************************************/

/*!\file
 * Implements operators on tImages.
 *
 * Usage:
 *   After doing
 *   \code
 *   include <gravis/tImage/operators.h>
 *   \endcode
 *   you can do pixelwise operations in images using overloaded functions.
 *
 *   For Example to rescale an image into [-1,1]:
 *   \code
 *   tImage< tRGB<double> > image;
 *   image.read('test.png');
 *   image -= 0.5;
 *   image *= 2.0;
 *   \endcode
 *
 */
#ifndef __GRAVIS__IMAGE_OPERATORS__
#define __GRAVIS__IMAGE_OPERATORS__

#include <src/jaz/gravis/tImage.h>

namespace gravis
{

  // \@{
  //! Inplace image - image operation
  template <class T1, class T2, class OP>
  void imageOpI(tImage<T1> &lhs, const tImage<T2> &rhs, const OP& op)
  {
    if ((lhs.cols() != rhs.cols()) ||
        (lhs.rows() != rhs.rows()))
      throw("Incompatible sizes.");
    T1* ldata = lhs.data();
    const T2* rdata = rhs.data();
    const T1* end = ldata + lhs.cols()*lhs.rows();
    for (; ldata<end; ++ldata, ++rdata)
      op(*ldata, *rdata);
  }

  //! Out of place image - image operation
  template <class T1, class T2, class OP>
  tImage<T1> imageOp(const tImage<T1> &rhs1, const tImage<T2> &rhs2, const OP& op)
  {
    tImage<T1> lhs(rhs2.cols(), rhs2.rows());
    T1* ldata = lhs.data();
    const T1* r1data = rhs1.data();
    const T2* r2data = rhs2.data();
    const T1* end = ldata + lhs.cols()*lhs.rows();
    for (; ldata<end; ++ldata, ++r1data, ++r2data)
      *ldata = op(*r1data, *r2data);
    return lhs;
  }

  //! Inplace image - scalar operation
  template <class T, class F, class OP>
  void imageOpI(tImage<T> &lhs, const F& rhs, const OP& op)
  {
    T* ldata = lhs.data();
    const T* end = ldata + lhs.cols()*lhs.rows();
    for (; ldata<end; ++ldata)
      op(*ldata, rhs);
  }

  //! Out of place image - scalar operation
  template <class T1, class F, class OP>
  tImage<T1> imageOp(const tImage<T1> &rhs1, const F& rhs2, const OP& op)
  {
    tImage<T1> lhs(rhs1.cols(), rhs1.rows());
    T1* ldata = lhs.data();
    const T1* rdata = rhs1.data();
    const T1* end = ldata + lhs.cols()*lhs.rows();
    for (; ldata<end; ++ldata, ++rdata)
      *ldata = op(*rdata, rhs2);
    return lhs;
  }

  //! Out of place scalar - image operation
  template <class T1, class F, class OP>
  tImage<T1> imageOp(const F& rhs1, const tImage<T1> &rhs2, const OP& op)
  {
    tImage<T1> lhs(rhs2.cols(), rhs2.rows());
    T1* ldata = lhs.data();
    const T1* rdata = rhs2.data();
    const T1* end = ldata + lhs.cols()*lhs.rows();
    for (; ldata<end; ++ldata, ++rdata)
      *ldata = op(rhs1, *rdata);
    return lhs;
  }

  //! Inplace unary operator application.
  template <class T, class OP>
  void imageOpI(tImage<T> &img, const OP& op)
  {
    T* data = img.data();
    const T* end = data + img.cols()*img.rows();
    for (; data<end; ++data)
      *data = op(*data);
  }

  //! Out of place unary operator application.
  template <class T, class OP>
  tImage<T> imageOp(const tImage<T> &img, const OP& op)
  {
    tImage<T> lhs = img.clone();
    imageOpI(lhs, op);
    return lhs;
  }

  // \@}

  namespace priv
  {
    //! Inplace Subtraction Functor
    template <class T1, class T2>  struct FunSubI
    {
      inline void operator()(T1& a, const T2& b) const
      {
        a -= b;
      }
    };
    //! Inplace Addition Functor
    template <class T1, class T2>  struct FunAddI
    {
      inline void operator()(T1& a, const T2& b) const
      {
        a += b;
      }
    };
    //! Inplace Multiplication Functor
    template <class T1, class T2>  struct FunMulI
    {
      inline void operator()(T1& a, const T2& b) const
      {
        a *= b;
      }
    };
    //! Inplace Division Functor
    template <class T1, class T2>  struct FunDivI
    {
      inline void operator()(T1& a, const T2& b) const
      {
        a /= b;
      }
    };

    //! Out of place Subtraction Functor
    template <class R, class T1, class T2>  struct FunSub
    {
      inline R operator()(const T1& a, const T2& b) const
      {
        return a - b;
      }
    };
    //! Out of place Addition Functor
    template <class R, class T1, class T2>  struct FunAdd
    {
      inline R operator()(const T1& a, const T2& b) const
      {
        return a + b;
      }
    };
    //! Out of place Multiplication Functor
    template <class R, class T1, class T2>  struct FunMul
    {
      inline R operator()(const T1& a, const T2& b) const
      {
        return a * b;
      }
    };
    //! Out of place Division Functor
    template <class R, class T1, class T2>  struct FunDiv
    {
      inline R operator()(const T1& a, const T2& b) const
      {
        return a / b;
      }
    };

    //! negation functor
    template <class T>  struct FunNeg
    {
      inline T operator()(const T& a) const
      {
        return -a;
      }
    };
    //! not functor
    template <class T>  struct FunNot
    {
      inline T operator()(const T& a) const
      {
        return !a;
      }
    };
    //! abs functor
    template <class T>  struct FunAbs
    {
      inline T operator()(const T& a) const
      {
        return abs(a);
      }
    };
  }

  //! Subtract one image from another inplace.
  template <class T1, class T2>
  void operator-=(tImage<T1> &lhs, const tImage<T2> &rhs)
  {
    priv::FunSubI<T1, T2> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Add one image to another inplace.
  template <class T1, class T2>
  void operator+=(tImage<T1> &lhs, const tImage<T2> &rhs)
  {
    priv::FunAddI<T1, T2> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Multiply an image with another inplace.
  template <class T1, class T2>
  void operator*=(tImage<T1> &lhs, const tImage<T2> &rhs)
  {
    priv::FunMulI<T1, T2> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Divide one image with another inplace.
  template <class T1, class T2>
  void operator/=(tImage<T1> &lhs, const tImage<T2> &rhs)
  {
    priv::FunDivI<T1, T2> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Subtract one image from another.
  template <class T1, class T2>
  tImage<T1> operator-(const tImage<T1> &rhs1, const tImage<T2> &rhs2)
  {
    priv::FunSub<T1, T1, T2> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Add one image to another
  template <class T1, class T2>
  tImage<T1> operator+(const tImage<T1> &rhs1, const tImage<T2> &rhs2)
  {
    priv::FunAdd<T1, T1, T2> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Multiply one image with another
  template <class T1, class T2>
  tImage<T1> operator*(const tImage<T1> &rhs1, const tImage<T2> &rhs2)
  {
    priv::FunMul<T1, T1, T2> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Divide one image by another
  template <class T1, class T2>
  tImage<T1> operator/(const tImage<T1> &rhs1, const tImage<T2> &rhs2)
  {
    priv::FunDiv<T1, T1, T2> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Subtract a scalar from an image inplace.
  template <class T, class F>
  void operator-=(tImage<T> &lhs, const F& rhs)
  {
    priv::FunSubI<T, F> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Add a scalar to an image inplace.
  template <class T, class F>
  void operator+=(tImage<T> &lhs, const F& rhs)
  {
    priv::FunAddI<T, F> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Multiply an image with a scalar inplace.
  template <class T, class F>
  void operator*=(tImage<T> &lhs, const F& rhs)
  {
    priv::FunMulI<T, F> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Divide an image by a scalar inplace.
  template <class T, class F>
  void operator/=(tImage<T> &lhs, const F& rhs)
  {
    priv::FunDivI<T, F> fun;
    imageOpI(lhs, rhs, fun);
  }

  //! Subtract an scalar from an image
  template <class T, class F>
  tImage<T> operator-(const tImage<T> &rhs1, const F& rhs2)
  {
    priv::FunSub<T, T, F> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Add a scalar to an image
  template <class T, class F>
  tImage<T> operator+(const tImage<T> &rhs1, const F& rhs2)
  {
    priv::FunAdd<T, T, F> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Multiply an image with a scalar
  template <class T, class F>
  tImage<T> operator*(const tImage<T> &rhs1, const F& rhs2)
  {
    priv::FunMul<T, T, F> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Divide an image by a scalar
  template <class T, class F>
  tImage<T> operator/(const tImage<T> &rhs1, const F& rhs2)
  {
    priv::FunDiv<T, T, F> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Left subtraction of scalar with image
  template <class T, class F>
  tImage<T> operator-(const F& rhs1, const tImage<T> &rhs2)
  {
    priv::FunSub<T, F, T> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Left addition of scalar to image
  template <class T, class F>
  tImage<T> operator+(const F& rhs1, const tImage<T> &rhs2)
  {
    priv::FunAdd<T, F, T> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Left multiply an image with a scalar
  template <class T, class F>
  tImage<T> operator*(const F& rhs1, const tImage<T> &rhs2)
  {
    priv::FunMul<T, F, T> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Left divide an image with a scalar
  template <class T, class F>
  tImage<T> operator/(const F& rhs1, const tImage<T> &rhs2)
  {
    priv::FunDiv<T, F, T> fun;
    return imageOp(rhs1, rhs2, fun);
  }

  //! Negate an image
  template <class T>
  tImage<T> operator-(const tImage<T> &img)
  {
    priv::FunNeg<T> fun;
    return imageOp(img, fun);
  }

  //! Calculate the absolute of an image
  template <class T>
  tImage<T> abs(const tImage<T> &img)
  {
    priv::FunAbs<T> fun;
    return imageOp(img, fun);
  }

  //! Calculate the absolute of an image inplace
  template <class T>
  void absI(tImage<T> &img)
  {
    priv::FunAbs<T> fun;
    return imageOpI(img, fun);
  }

  //! Negate an image
  template <class T>
  tImage<T> operator!(const tImage<T> &img)
  {
    priv::FunNot<T> fun;
    return imageOp(img, fun);
  }

}

#endif
