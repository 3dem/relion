/******************************************************************************
**        Title: Interpolation for tImage
**  Description: Interpolated Image Access
**
**       Author: Brian Amberg, 2006
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#ifndef __GRAVIS__IMAGE_INTERPOLATION__
#define __GRAVIS__IMAGE_INTERPOLATION__

#include "../tImage.h"
#include "../tImage/traits.h"
#include <cmath>

namespace gravis
{
  namespace interpolation
  {

    /**
     * Nearest Neighbour Image Access
     *
     * Usage:
     *   image.interpolate<interpolation::NearestNeighbour>(20.2, 20.4)
     **/
    struct NearestNeighbour
    {
      private:
        template <class T> inline static int round(const T& v)
        {
          return v < 0 ? int(v-T(0.5)) : int(v+T(0.5));
        }

      public:
        template <class T, class F>
        static inline
        T getPixel(const tImage<T>& image, const F& x, const F& y)
        {
          const int x_i = round(x);
          const int y_i = round(y);
          return image.access(x_i, y_i);
        }
    };

    /**
     * Linearly interpolated image access
     *
     * Usage:
     *   image.interpolate<interpolation::Linear>(20.2, 20.4)
     **/
    struct Linear
    {
      private:
        template <class T> inline static int round(const T& v)
        {
          return v < 0 ? int(v-T(0.5)) : int(v+T(0.5));
        }

      public:
        template <class T, class F>
        static inline
        T getPixel(const tImage<T>& image, const F& x, const F& y)
        {
          const int x_i = int(floor(x));
          const int y_i = int(floor(y));
          const F   dx  = x-x_i;
          const F   dy  = y-y_i;
          return
            (image.access(x_i   , y_i  ) * (F(1.0)-dx) +
             image.access(x_i+1 , y_i  ) * (       dx)) * (F(1.0)-dy) +
            (image.access(x_i   , y_i+1) * (F(1.0)-dx) +
             image.access(x_i+1 , y_i+1) * (       dx)) * (       dy);
        }
    };


    /**
     * Cubic interpolated image access
     *
     * Usage:
     *   image.interpolate<interpolation::Cubic>(20.2, 20.4)
     **/
    struct Cubic
    {
      private:
        template <class T> inline static int round(const T& v)
        {
          return v < 0 ? int(v-T(0.5)) : int(v+T(0.5));
        }

        template <class T, class F>
        static inline
        T cubicInterpolation(const T& a, const T& b, const T& c, const T& d, const F& x)
        {
          const T p = (d - c) - (a - b);
          const T q = (a - b) - p;
          const T r = c - a;

          return p*(x*x*x) + q*(x*x) + r*x + b;
        }

      public:
        template <class T, class F>
        static inline
        T getPixel(const tImage<T>& image, const F& x, const F& y)
        {
          int x_i = int(floor(x));
          int y_i = int(floor(y));
          const F dx = x-x_i;
          const F dy = y-y_i;
          return cubicInterpolation(
                   cubicInterpolation(
                     image.access(x_i-1, y_i-1),
                     image.access(x_i,   y_i-1),
                     image.access(x_i+1, y_i-1),
                     image.access(x_i+2, y_i-1), dx),
                   cubicInterpolation(
                     image.access(x_i-1, y_i  ),
                     image.access(x_i,   y_i  ),
                     image.access(x_i+1, y_i  ),
                     image.access(x_i+2, y_i  ), dx),
                   cubicInterpolation(
                     image.access(x_i-1, y_i+1),
                     image.access(x_i,   y_i+1),
                     image.access(x_i+1, y_i+1),
                     image.access(x_i+2, y_i+1), dx),
                   cubicInterpolation(
                     image.access(x_i-1, y_i+2),
                     image.access(x_i,   y_i+2),
                     image.access(x_i+1, y_i+2),
                     image.access(x_i+2, y_i+2), dx), dy);
        }
    };

    /**
     * Wrapper around another interpolation method, that first scales textures from [0,1]x[0,1] to [0,width-1]x[0,height-1]
     *
     * Usage:
     *   image.interpolate<interpolation::NearestNeighbourTextureCoordinate>(0.2, 0.4)
     **/
    template <class AccessMethod>
    struct TextureCoordinateAccess
    {
      template <class T, class F>
      static inline
      T getPixel(const tImage<T>& image, const F& x, const F& y)
      {
        return AccessMethod::getPixel(image, x*F(image.cols()-1), y*F(image.rows()-1));
      }
    };

    /**
     * Nearest Neighbour interpolated access to the image, using coordinates in [0,1]
     *
     * Usage:
     *   image.interpolate<interpolation::NearestNeighbourTextureCoordinate>(0.2, 0.4)
     **/
    typedef TextureCoordinateAccess<NearestNeighbour> NearestNeighbourTextureCoordinate;

    /**
     * Linear interpolated access to the image, using coordinates in [0,1]
     *
     * Usage:
     *   image.interpolate<interpolation::NearestNeighbourTextureCoordinate>(0.2, 0.4)
     **/
    typedef TextureCoordinateAccess<Linear> LinearTextureCoordinate;

    /**
     * Cubic interpolated access to the image, using coordinates in [0,1]
     *
     * Usage:
     *   image.interpolate<interpolation::NearestNeighbourTextureCoordinate>(0.2, 0.4)
     **/
    typedef TextureCoordinateAccess<Cubic> CubicTextureCoordinate;

  }
}


#endif
