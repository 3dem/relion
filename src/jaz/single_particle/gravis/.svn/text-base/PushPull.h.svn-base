/*============================================================================*/
/**
 *         @file PushPull.h
 *
 *        @brief Push-Pull interpolation class
 *
 *         @date 1 Aug 2012
 *      @authors Jasenko Zivanov\n
 *               jasenko.zivanov@unibas.ch\n
 *               University of Basel, Switzerland
 */
/*============================================================================*/

#ifndef PUSH_PULL_H
#define PUSH_PULL_H

#include <vector>
#include "tImage.h"

namespace gravis
{
  /** Push-Pull interpolation, interpolation on the image pyramid filling holes in an image */
  class PushPull
  {
    public:
      /** Execute a push-pull interpolation (image pyramid hole filler) on the image
        * \param img Image to interpolate, do not use an alpha channel, use the mask to specify holes
        * \param mask Mask indicating which regions are holes (0: fill, 1: keep - as in the alpha channel)
       */
      template <class T>
      static gravis::tImage<T> interpolate(
          gravis::tImage<T> img,
          gravis::tImage<float> mask,
          int minSize = 1);

    private:
      template <class T>
      static gravis::tImage<T> shrink2(gravis::tImage<T> img);

      template <class T>
      static gravis::tImage<T> shrink2(gravis::tImage<T> img, gravis::tImage<float> mask);

      template <class T>
      static gravis::tImage<T> grow2(gravis::tImage<T> img);

      template <class T>
      static gravis::tImage<T> blur3x3(gravis::tImage<T> img);
  };

  template <class T>
  gravis::tImage<T> PushPull :: interpolate(
      gravis::tImage<T> img,
      gravis::tImage<float> mask,
      int minSize)
  {
    const int w = img.cols();
    const int h = img.rows();

    gravis::tImage<T> out(w,h);

    std::vector<gravis::tImage<T> > pyramid(0);
    std::vector<gravis::tImage<float> > maskPyramid(0);

    pyramid.push_back(img);
    maskPyramid.push_back(mask);

    gravis::tImage<T> ci = img;
    gravis::tImage<float> cm = mask;

    while ((int)ci.rows() > minSize && (int)ci.cols() > minSize)
    {
      ci = shrink2(ci,cm);
      cm = shrink2(cm);

      pyramid.push_back(ci);
      maskPyramid.push_back(cm);
    }

    maskPyramid[pyramid.size()-2].fill(1.f);

    for (int i = (int)pyramid.size() - 2; i >= 0; i--)
    {
      gravis::tImage<T> pi0 = pyramid[i];
      gravis::tImage<T> pi1 = grow2(pyramid[i+1]);
      gravis::tImage<float> pm = maskPyramid[i];

      for (int y = 0; y < (int)pi0.rows(); y++)
        for (int x = 0; x < (int)pi0.cols(); x++)
        {
          pi0(x,y) = pm(x,y)*pi0(x,y) + (1.f - pm(x,y))*pi1(x,y);
        }
    }

    return pyramid[0];
  }

  template <class T>
  gravis::tImage<T> PushPull :: shrink2(gravis::tImage<T> img, gravis::tImage<float> mask)
  {
    const int w = img.cols();
    const int h = img.rows();

    const int w2 = (int)(ceil(w/2.f));
    const int h2 = (int)(ceil(h/2.f));

    gravis::tImage<T> out(w2,h2);
    gravis::tImage<float> weight(w2,h2);

    out.fill(T(0.f));
    weight.fill(0.f);

    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++)
      {
        out(x/2, y/2) += mask(x,y) * img(x,y);
        weight(x/2, y/2) += mask(x,y);
      }

    for (int i = 0; i < (w/2)*(h/2); i++)
    {
      if (weight[i] > 0.f) out[i] /= weight[i];
    }

    return out;
  }

  template <class T>
  gravis::tImage<T> PushPull :: shrink2(gravis::tImage<T> img)
  {
    const int w = img.cols();
    const int h = img.rows();

    const int w2 = (int)(ceil(w/2.f));
    const int h2 = (int)(ceil(h/2.f));

    gravis::tImage<T> out(w2,h2);
    gravis::tImage<float> weight(w2,h2);

    out.fill(T(0.f));
    weight.fill(0.f);

    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++)
      {
        out(x/2, y/2) += img(x,y);
        weight(x/2, y/2)++;
      }

    for (int i = 0; i < w2*h2; i++)
    {
      if (weight[i] > 0.f) out[i] /= weight[i];
    }

    return out;
  }

  template <class T>
  gravis::tImage<T> PushPull :: grow2(gravis::tImage<T> img)
  {
    const int w = img.cols();
    const int h = img.rows();

    gravis::tImage<T> out(2*w,2*h);

    for (int y = 0; y < 2*h; y++)
      for (int x = 0; x < 2*w; x++)
      {
        out(x, y) = img(x/2,y/2);
      }

    return blur3x3(out);
  }

  template <class T>
  gravis::tImage<T> PushPull :: blur3x3(gravis::tImage<T> img)
  {
    const int w = img.cols();
    const int h = img.rows();

    gravis::tImage<T> out(w,h);
    gravis::tImage<float> weight(w,h);

    out.fill(T(0.f));
    weight.fill(0.f);

    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++)
      {
        if (x > 0 && y > 0)
        {
          out(x,y) += img(x-1,y-1)/4.f;
          weight(x,y) += 1/4.f;
        }
        if (x > 0)
        {
          out(x,y) += img(x-1,y)/2.f;
          weight(x,y) += 1/2.f;
        }
        if (x > 0 && y < h-1)
        {
          out(x,y) += img(x-1,y+1)/4.f;
          weight(x,y) += 1/4.f;
        }

        if (y > 0)
        {
          out(x,y) += img(x,y-1)/2.f;
          weight(x,y) += 1/2.f;
        }
        {
          out(x,y) += img(x,y);
          weight(x,y)++;
        }
        if (y < h-1)
        {
          out(x,y) += img(x,y+1)/2.f;
          weight(x,y) += 1/2.f;
        }

        if (x < w-1 && y > 0)
        {
          out(x,y) += img(x+1,y-1)/4.f;
          weight(x,y) += 1/4.f;
        }
        if (x < w-1)
        {
          out(x,y) += img(x+1,y)/2.f;
          weight(x,y) += 1/2.f;
        }
        if (x < w-1 && y < h-1)
        {
          out(x,y) += img(x+1,y+1)/4.f;
          weight(x,y) += 1/4.f;
        }
      }

    for (int i = 0; i < w*h; i++)
    {
      out[i] /= weight[i];
    }

    return out;
  }
}

#endif
