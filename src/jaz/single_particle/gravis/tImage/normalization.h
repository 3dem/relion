/*************************************************************************//*!
*        Title: tImageNormalization.h
*  Description: Implements image normalization
*
*       Author: Brian Schroeder, 2006
*               Computer Science Department, University Basel (CH)
****************************************************************************/

/*!\file
 * Implements image normalization on scalar, rgb, rgba and vector 1,2,3 entries.
 *
 * Usage:
 *   After doing
 *   \code
 *   include <gravis/tImage/normalization.h>
 *   \endcode
 *   you can use normalize() and normalizeI() on the images.
 *
 */
#ifndef __GRAVIS__TIMAGE_NORMALIZATION__
#define __GRAVIS__TIMAGE_NORMALIZATION__

#include <src/jaz/gravis/tImage.h>
#include <src/jaz/gravis/tImage/traits.h>
#include <src/jaz/gravis/tRGB.h>
#include <src/jaz/gravis/tRGBA.h>
#include <src/jaz/gravis/tRGB_A.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t4Vector.h>

namespace gravis
{

  //! @cond INTERN
  namespace priv
  {
    //! Scalar minmax
    template <class T>
    struct FunMinMaxS
    {
      void operator()(T& min, T& max, const T& p) const
      {
        min = p;
        max = p;
      };
    };

    //! access to .r .g .b
    template <class T, class RGB>
    struct FunMinMaxRGB
    {
      void operator()(T& min, T& max, const RGB& p) const
      {
        min = std::min(std::min(p.r, p.g), p.b);
        max = std::max(std::max(p.r, p.g), p.b);
      };
    };

    //! access to .g
    template <class T, class GA>
    struct FunMinMaxGA
    {
      void operator()(T& min, T& max, const GA& p) const
      {
        min = p.g;
        max = p.g;
      };
    };

    //! access to .r .g .b
    template <class RGB>
    struct FunMinMaxRGBChannel
    {
      void operator()(RGB& min, RGB& max, const RGB& p) const
      {
        min.r = std::min(min.r,p.r);
        min.g = std::min(min.g,p.g);
        min.b = std::min(min.b,p.b);

        max.r = std::max(max.r,p.r);
        max.g = std::max(max.g,p.g);
        max.b = std::max(max.b,p.b);
      };
    };

    //! Things accessible by operator[]
    template <class T, class V, int i>
    struct FunMinMaxV
    {
      void operator()(T& min, T& max, const V& p) const
      {
        min = p[0];
        max = p[0];
        for (int j=1; j<i; j++)
        {
          min = std::min(min, p[j]);
          max = std::max(max, p[j]);
        };
      }
    };
  };

  //! Inplace image normalization
  template <class T, class FUN>
  void normalizeI(tImage< T > &in, const FUN& fun)
  {

    typedef typename tImageTraits< T >::Scalar_t S;
    typedef typename tImageTraits< T >::Float_t  F;

    T* data = in.data();
    if(!data) return;
    const T* end = data + in.rows() * in.cols();
    S min, max;
    fun(min, max, *data);
    for (const T* p = data; p<end; p++)
    {
      S _min, _max;
      fun(_min, _max, *p);
      min = std::min(min, _min);
      max = std::max(max, _max);
    }
    const F c = F(1) / (F(max) - F(min));
    const S o = -min;
    for (T* p = data; p<end; p++)
    {
      *p += o;
      *p = *p*c;
    }
  };

  //! Out of place image normalization
  template <class T, class FUN>
  tImage< T > normalize(const tImage< T > &in, const FUN& fun)
  {
    tImage< T > lhs = in.clone();
    normalizeI<T, FUN>(lhs, fun);
    return lhs;
  };
  //! @endcond

  //! Scalar Image inplace normalization
  // Scales and offsets the pixel values, such that the max and min over all pixels is 0 and 1 respectively.
  template <class T> void normalizeI(tImage< T > &in)
  {
    priv::FunMinMaxS<T> fun;
    normalizeI< T, priv::FunMinMaxS<T> >(in, fun);
  }
  //! Scalar Image normalization
  // Scales and offsets the pixel values, such that the max and min over all pixels is 0 and 1 respectively.
  template <class T> tImage< T > normalize(const tImage< T > &in)
  {
    priv::FunMinMaxS<T> fun;
    return normalize< T, priv::FunMinMaxS<T> >(in, fun);
  }

  //! Gray_A Image inplace normalization
  // Scales and offsets the pixel values, such that the max and min over all pixels is 0 and 1 respectively.
  template <class T> void normalizeI(tImage<tGray_A<T> > &in)
  {
    priv::FunMinMaxGA<T, tGray_A<T> > fun;
    normalizeI< tGray_A<T>, priv::FunMinMaxGA<T, tGray_A<T> > >(in, fun);
  }
  //! Gray_A Image normalization
  // Scales and offsets the pixel values, such that the max and min over all pixels is 0 and 1 respectively.
  template <class T> tImage< tGray_A<T> > normalize(const tImage< tGray_A<T> > &in)
  {
    priv::FunMinMaxGA<T, tGray_A<T> > fun;
    return normalize< T, priv::FunMinMaxGA<T, tGray_A<T> > >(in, fun);
  }

  //! RGB Image inplace normalization
  // Scales and offsets all channels simultaneously, such that the max and min over all pixels and channels is 0 and 1 respectively.
  template <class T> void normalizeI(tImage< tRGB< T > > &in)
  {
    priv::FunMinMaxRGB< T, tRGB<T> > fun;
    normalizeI< tRGB< T >, priv::FunMinMaxRGB< T, tRGB<T> > >(in, fun);
  }
  //! RGB Image normalization
  // Scales and offsets all channels simultaneously, such that the max and min over all pixels and channels is 0 and 1 respectively.
  template <class T> tImage< tRGB< T > > normalize(const tImage< tRGB< T > > &in)
  {
    priv::FunMinMaxRGB< T, tRGB<T> > fun;
    return normalize< tRGB< T >, priv::FunMinMaxRGB< T, tRGB<T> > >(in, fun);
  }

  //! RGBA Image inplace normalization
  // Scales and offsets all channels simultaneously, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> void normalizeI(tImage< tRGBA< T > > &in)
  {
    priv::FunMinMaxRGB< T, tRGBA<T> > fun;
    normalizeI< tRGBA< T >, priv::FunMinMaxRGB< T, tRGBA<T> > >(in, fun);
  }
  //! RGBA Image normalization
  // Scales and offsets all channels simultaneously, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> tImage< tRGBA< T > > normalize(tImage< tRGBA< T > > &in)
  {
    priv::FunMinMaxRGB< T, tRGBA<T> > fun;
    return normalize< tRGBA< T >, priv::FunMinMaxRGB< T, tRGBA<T> > >(in, fun);
  }
  //! RGB_A Image inplace normalization
  // Scales and offsets all channels simultaneously, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> void normalizeI(tImage< tRGB_A< T > > &in)
  {
    priv::FunMinMaxRGB< T, tRGB_A<T> > fun;
    normalizeI< tRGB_A< T >, priv::FunMinMaxRGB< T, tRGB_A<T> > >(in, fun);
  }
  //! RGB_A Image normalization
  // Scales and offsets all channels simultaneously, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> tImage< tRGB_A< T > > normalize(tImage< tRGB_A< T > > &in)
  {
    priv::FunMinMaxRGB< T, tRGB_A<T> > fun;
    return normalize< tRGB_A< T >, priv::FunMinMaxRGB< T, tRGB_A<T> > >(in, fun);
  }


  //! RGB Image inplace normalization
  // Scales and offsets all channels seperatly, such that the max and min over all pixels and channels is 0 and 1 respectively.
  template <class T> void normalizeIC(tImage< tRGB< T > > &in)
  {
    priv::FunMinMaxRGBChannel<tRGB<T> > fun;
    normalizeI< tRGB< T >, priv::FunMinMaxRGBChannel<tRGB<T> > >(in, fun);
  }
  //! RGB Image normalization
  // Scales and offsets all channels seperatly, such that the max and min over all pixels and channels is 0 and 1 respectively.
  template <class T> tImage< tRGB< T > > normalizeC(const tImage< tRGB< T > > &in)
  {
    priv::FunMinMaxRGBChannel<tRGB<T> > fun;
    return normalize< tRGB< T >, priv::FunMinMaxRGBChannel<tRGB<T> > >(in, fun);
  }

  //! RGBA Image inplace normalization
  // Scales and offsets all channels seperatly, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> void normalizeIC(tImage< tRGBA< T > > &in)
  {
    priv::FunMinMaxRGBChannel<tRGBA<T> > fun;
    normalizeI< tRGBA< T >, priv::FunMinMaxRGBChannel<tRGBA<T> > >(in, fun);
  }
  //! RGBA Image normalization
  // Scales and offsets all channels seperatly, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> tImage< tRGBA< T > > normalizeC(tImage< tRGBA< T > > &in)
  {
    priv::FunMinMaxRGBChannel<tRGBA<T> > fun;
    return normalize< tRGBA< T >, priv::FunMinMaxRGBChannel<tRGBA<T> > >(in, fun);
  }

  //! RGB_A Image inplace normalization
  // Scales and offsets all channels seperatly, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> void normalizeIC(tImage< tRGB_A< T > > &in)
  {
    priv::FunMinMaxRGBChannel<tRGB_A<T> > fun;
    normalizeI< tRGB_A< T >, priv::FunMinMaxRGBChannel<tRGB_A<T> > >(in, fun);
  }
  //! RGB_A Image normalization
  // Scales and offsets all channels seperatly, such that the max and min over all pixels and channels is 0 and 1 respectively. Ignores the alpha channel.
  template <class T> tImage< tRGB_A< T > > normalizeC(tImage< tRGB_A< T > > &in)
  {
    priv::FunMinMaxRGBChannel<tRGB_A<T> > fun;
    return normalize< tRGB_A< T >, priv::FunMinMaxRGBChannel<tRGB_A<T> > >(in, fun);
  }

  //! t2Vector Image inplace normalization
  // Scales and offsets all dimensions simultaneously, such that the max and min over all pixels and dimensions is 0 and 1 respectively.
  template <class T> void normalizeI(tImage< t2Vector< T > > &in)
  {
    priv::FunMinMaxV< T, t2Vector<T>, 2> fun;
    normalizeI< t2Vector< T >, priv::FunMinMaxV< T, t2Vector<T>, 2> >(in, fun);
  }
  //! t 2Vector Image normalization
  // Scales and offsets all dimensions simultaneously, such that the max and min over all pixels and dimensions is 0 and 1 respectively.
  template <class T> tImage< t2Vector< T > > normalize(const tImage< t2Vector< T > > &in)
  {
    priv::FunMinMaxV< T, t2Vector<T>, 2> fun;
    return normalize< t2Vector< T >, priv::FunMinMaxV< T, t2Vector<T>, 2> >(in, fun);
  }

  //! t3Vector Image inplace normalization
  // Scales and offsets all dimensions simultaneously, such that the max and min over all pixels and dimensions is 0 and 1 respectively.
  template <class T> void normalizeI(tImage< t3Vector< T > > &in)
  {
    priv::FunMinMaxV< T, t3Vector<T>, 2> fun;
    normalizeI< t3Vector< T >, priv::FunMinMaxV< T, t3Vector<T>, 2> >(in, fun);
  }
  //! t3Vector Image normalization
  // Scales and offsets all dimensions simultaneously, such that the max and min over all pixels and dimensions is 0 and 1 respectively.
  template <class T> tImage< t3Vector< T > > normalize(const tImage< t3Vector< T > > &in)
  {
    priv::FunMinMaxV< T, t3Vector<T>, 3> fun;
    return normalize< t3Vector< T >, priv::FunMinMaxV< T, t3Vector<T>, 3> >(in, fun);
  }

  //! t4Vector Image inplace normalization
  // Scales and offsets all dimensions simultaneously, such that the max and min over all pixels and dimensions is 0 and 1 respectively.
  template <class T> void normalizeI(tImage< t4Vector< T > > &in)
  {
    priv::FunMinMaxV< T, t4Vector<T>, 2> fun;
    normalizeI< t4Vector< T >, priv::FunMinMaxV< T, t4Vector<T>, 2> >(in, fun);
  }
  //! t4Vector Image normalization
  // Scales and offsets all dimensions simultaneously, such that the max and min over all pixels and dimensions is 0 and 1 respectively.
  template <class T> tImage< t4Vector< T > > normalize(const tImage< t4Vector< T > > &in)
  {
    priv::FunMinMaxV< T, t4Vector<T>, 4> fun;
    return normalize< t4Vector< T >, priv::FunMinMaxV< T, t4Vector<T>, 4> >(in, fun);
  }
};

#endif
