/*************************************************************************//*!
*        Title: tImage/distance_transform.h
*  Description: Implements drawing operations on images without the need
*               for the libRender library.
*
*       Author: Brian Amberg, 2006-2007
*               Computer Science Department, University Basel (CH)
****************************************************************************/
#ifndef __GRAVIS__TIMAGE__DISTANCE_TRANSFORM__
#define __GRAVIS__TIMAGE__DISTANCE_TRANSFORM__

#include "../tImage/traits.h"

namespace gravis
{

  template <class T> inline static T sqr(const T& a)
  {
    return a*a;
  }

  /**
   * Binarize an image using an intensity threshold
   **/
  template <class T, class F>
  static inline
  void intensityThreshold(tImage<bool> &out, const tImage<T> &in, const F& threshold = 0)
  {
    out.resize(in.cols(), in.rows());
    for (size_t i=0; i<in.size(); ++i) out[i] = in[i].intensity() > threshold;
  };

  /**
   * Binarize an image using a threshold on the alpha channel
   **/
  template <class T, class F>
  static inline
  void alphaThreshold(tImage<bool> &out, const tImage<T> &in, const F& threshold = typename tImageTraits<T>::Float_t(0))
  {
    out.resize(in.cols(), in.rows());
    for (size_t i=0; i<in.size(); ++i) out[i] = in[i].a > threshold;
  };

  /**
   * w=1, iterations=2 gives perfect results, while still being fast
   **/
  template <class F, class Pixel>
  inline static
  void distanceTransform(tImage< F > &out, const tImage< Pixel > &in, const int& w=3, const size_t& iterations=3)
  {
    const int W=in.cols();
    const int H=in.rows();

    out.resize(W,H);

    // Temporary Memory (parents)
    tImage< t2Vector< int > > p(W, H);

    tImage<F> dm(2*w+1, 2*w+1); // distances
    for (int x=0; x<2*w+1; ++x)
      for (int y=0; y<2*w+1; ++y)
        dm(x,y) = sqrt(sqr(y-w)+sqr(x-w));

    // initialize
    out.fill( in.cols()+in.rows() );
    p.fill( t2Vector<int>(0,0) );

    //initialize immediate interior elements
    for (int x=0; x<W; x++)
    {
      for (int y=0; y<H; y++)
      {
        /*
        if ( in(x,y) ) {
          out(x,y) = 0; // = 0 (distance trafo)
          p(x,y)[0] = x;
          p(x,y)[1] = y;
        }*/

        out(x,y) = in(x,y); // = 0 (distance trafo)
        p(x,y)[0] = x;
        p(x,y)[1] = y;
      }
    }

    for (size_t iteration=0; iteration<iterations; ++iteration)
    {

#define distanceTransformApplyRegion                                                                       \
        for (int l=0; l<2*w+1; ++l)                                                                        \
        for (int k=0; k<2*w+1; ++k) {                                                                      \
          if (out(x+l-w,y+k-w)+dm(l,k) < out(x,y)) {                                                       \
            out(x,y) = sqrt(sqr(p(x+l-w,y+k-w)[1] - y) + sqr(p(x+l-w,y+k-w)[0] - x));                      \
            p(x,y) = p(x+l-w,y+k-w); }                                                                     \
        }
#define distanceTransformApplyRegionSave                                                                   \
        for (int l=0; l<2*w+1; ++l)                                                                        \
        for (int k=0; k<2*w+1; ++k) {                                                                      \
          if (out.access(x+l-w,y+k-w)+dm(l,k) < out(x,y)) {                                                \
            out(x,y) = sqrt(sqr(p.access(x+l-w,y+k-w)[1] - y) + sqr(p.access(x+l-w,y+k-w)[0] - x));        \
            p(x,y) = p.access(x+l-w,y+k-w);                                                                \
          }                                                                                                \
        }
      //perform the forward pass
      for (int y=w; y<H-w; y++) for (int x=w; x<W-w; x++) distanceTransformApplyRegion;

      //Make the border
      for (int y=w-1; y>-1; y--) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=H-w; y<H;  y++) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=w-1; x>-1; x--) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=W-w; x<W;  x++) distanceTransformApplyRegionSave;

      // backward pass
      for (int y=H-w-1; y>=w; y--) for (int x=W-w-1; x>=w; x--) distanceTransformApplyRegion;

      //Make the border
      for (int y=w-1; y>-1; y--) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=H-w; y<H;  y++) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=w-1; x>-1; x--) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=W-w; x<W;  x++) distanceTransformApplyRegionSave;

#undef distanceTransformApplyRegion
#undef distanceTransformApplyRegionSave
    }

  }

  /**
   * w=1, iterations=2 gives perfect results, while still being fast
   **/
  template <class F, class Pixel>
  inline static
  void distanceTransformSq(tImage< F > &out, const tImage< Pixel > &in, const int& w=3, const size_t& iterations=3)
  {
    const int W=in.cols();
    const int H=in.rows();

    out.resize(W,H);

    // Temporary Memory (parents)
    tImage< t2Vector< int > > p(W, H);

    tImage<F> dm(2*w+1, 2*w+1); // distances
    for (int x=0; x<2*w+1; ++x)
      for (int y=0; y<2*w+1; ++y)
        dm(x,y) = sqr(y-w)+sqr(x-w);

    // initialize
    out.fill( in.cols()+in.rows() );
    p.fill( t2Vector<int>(0,0) );

    //initialize immediate interior elements
    for (int x=0; x<W; x++)
    {
      for (int y=0; y<H; y++)
      {
        /*
        if ( in(x,y) ) {
          out(x,y) = 0; // = 0 (distance trafo)
          p(x,y)[0] = x;
          p(x,y)[1] = y;
        }*/

        out(x,y) = in(x,y); // = 0 (distance trafo)
        p(x,y)[0] = x;
        p(x,y)[1] = y;
      }
    }

    for (size_t iteration=0; iteration<iterations; ++iteration)
    {

#define distanceTransformApplyRegion                                                                       \
        for (int l=0; l<2*w+1; ++l)                                                                        \
        for (int k=0; k<2*w+1; ++k) {                                                                      \
          if (out(x+l-w,y+k-w)+dm(l,k) < out(x,y)) {                                                       \
            out(x,y) = sqrt(sqr(p(x+l-w,y+k-w)[1] - y) + sqr(p(x+l-w,y+k-w)[0] - x));                      \
            p(x,y) = p(x+l-w,y+k-w); }                                                                     \
        }
#define distanceTransformApplyRegionSave                                                                   \
        for (int l=0; l<2*w+1; ++l)                                                                        \
        for (int k=0; k<2*w+1; ++k) {                                                                      \
          if (out.access(x+l-w,y+k-w)+dm(l,k) < out(x,y)) {                                                \
            out(x,y) = sqrt(sqr(p.access(x+l-w,y+k-w)[1] - y) + sqr(p.access(x+l-w,y+k-w)[0] - x));        \
            p(x,y) = p.access(x+l-w,y+k-w);                                                                \
          }                                                                                                \
        }
      //perform the forward pass
      for (int y=w; y<H-w; y++) for (int x=w; x<W-w; x++) distanceTransformApplyRegion;

      //Make the border
      for (int y=w-1; y>-1; y--) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=H-w; y<H;  y++) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=w-1; x>-1; x--) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=W-w; x<W;  x++) distanceTransformApplyRegionSave;

      // backward pass
      for (int y=H-w-1; y>=w; y--) for (int x=W-w-1; x>=w; x--) distanceTransformApplyRegion;

      //Make the border
      for (int y=w-1; y>-1; y--) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=H-w; y<H;  y++) for (int x=0;   x<W;  x++) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=w-1; x>-1; x--) distanceTransformApplyRegionSave;
      for (int y=0;   y<H;  y++) for (int x=W-w; x<W;  x++) distanceTransformApplyRegionSave;

#undef distanceTransformApplyRegion
#undef distanceTransformApplyRegionSave
    }

  }
}

#endif
