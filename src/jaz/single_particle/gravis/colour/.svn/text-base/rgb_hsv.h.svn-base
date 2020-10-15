#ifndef __GRAVIS__COLOUR__HSV__
#define __GRAVIS__COLOUR__HSV__

#include "../tRGBA.h"

// Author: Reinhard Knothe
// Extracted: Brian Schroeder

/**********************************************************************/
/* hsv2rgb:                                                           */
/**********************************************************************/
/* <-- rgb coordinates of hsv color, rgb in [0,1]^3                   */
/* --> hsv color coordinates h in [0, 360), s in [0,1], v in [0,1]    */
/**********************************************************************/
/* Converts hsv to rgb                                                */
/**********************************************************************/

template <class T>
gravis::tRGBA<T> hsv2rgb(T h, T s, T v, T a = T(1))
{

  if (s==0 && h==-1) return gravis::tRGBA<T>(v, v, v, a);
  else
  {
    if (h==360.0) h = 0.0;
    h /= 60.f;

    int i = (int)floor(h);

    T f = h - (T)i;

    T m = v * (1. - s);
    T n = v * (1. - s*f);
    T k = v * (1. - s*(1. - f));

    if (i==0) return gravis::tRGBA<T>(v, k, m, a);
    if (i==1) return gravis::tRGBA<T>(n, v, m, a);
    if (i==2) return gravis::tRGBA<T>(m, v, k, a);
    if (i==3) return gravis::tRGBA<T>(m, n, v, a);
    if (i==4) return gravis::tRGBA<T>(k, m, v, a);
    if (i==5) return gravis::tRGBA<T>(v, m, n, a);
  }

  return gravis::tRGBA<T>(v);
}

template <class T>
gravis::tRGBA<T> hsvScale(T dist, T min = 0.0, T max = 1.0, T a = 1.0)
{

  T h = 240.0 - (300.0 * (dist-min)/(max-min));
  if (h>360.0) h = 360.0;
  if (h<0.0) h = 0.0;

  return gravis::tRGBA<T>(hsv2rgba(h,(T) 1,(T) 1,a));

}

template <class T>
gravis::tRGBA<T> hsvBlueToRed(T dist, T min = 0.0, T max = 1.0, T a = 1.0)
{

  T d = (dist-min)/(max-min);
  if (d<0.0) d = 0.0;
  if (d>1.0) d = 1.0;
  T h = 240 - (240.0 * d);

  return gravis::tRGBA<T>(hsv2rgba(h,T(1), T(1), a));

}
template <class T>
gravis::tRGBA<T> hsvGreenToRed(T dist, T min = 0.0, T max = 1.0, T a = 1.0)
{

  T d = (dist-min)/(max-min);
  if (d<0.0) d = 0.0;
  if (d>1.0) d = 1.0;
  T h = (240.0 * d + 120);

  return gravis::tRGBA<T>(hsv2rgba(h,T(1), T(1), a));

}

#endif
