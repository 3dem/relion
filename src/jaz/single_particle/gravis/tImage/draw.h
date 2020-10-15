/*************************************************************************//*!
*        Title: tImage/draw.h
*  Description: Implements drawing operations on images without the need
*               for the libRender library.
*
*       Author: Brian Amberg, 2006-2007
*               Computer Science Department, University Basel (CH)
****************************************************************************/
#ifndef __GRAVIS__TIMAGE__DRAW__
#define __GRAVIS__TIMAGE__DRAW__

#include "../tImage.h"
#include "../tImage/traits.h"
#include "../t2Vector.h"
#include "../t3Vector.h"

namespace gravis
{

  namespace image
  {

    /**
     * Draws a line on an image, interpolating between the colours at the endpoint
     **/
    template <class F, class Pixel>
    static inline
    void draw_line(gravis::tImage<Pixel> &I,
                   const gravis::t2Vector<F> &a,  const gravis::t2Vector<F> &b,
                   const Pixel&               c1, const Pixel&               c2)
    {

      typedef typename tImageTraits<Pixel>::Float_t      ImageFloat;

      // Bounding Box
      const F min_x_f = std::max<F>(F(0),             std::min(a[0], b[0]));
      const F max_x_f = std::min<F>(F(I.cols()-1),    std::max(a[0], b[0]));
      const F min_y_f = std::max<F>(F(0),             std::min(a[1], b[1]));
      const F max_y_f = std::min<F>(F(I.rows()-1),    std::max(a[1], b[1]));

      // Bounding Box
      const size_t min_x = int( ceil(  min_x_f ) );
      const size_t max_x = int( floor( max_x_f ) );
      const size_t min_y = int( ceil(  min_y_f ) );
      const size_t max_y = int( floor( max_y_f ) );

      F alpha;
      // Loop over bounding box
      if (max_x - min_x > max_y - min_y)
      {

        if (a.x < b.x)
        {
          for (size_t x=min_x; x<=max_x; ++x)
          {
            alpha =  (F(x) - a.x) / (b.x - a.x);
            const size_t y = size_t(a.y + alpha * (b.y - a.y) + F(0.5));
            I(x, y) = ImageFloat(F(1) - alpha) * c1 + ImageFloat(alpha) * c2;
          }
        }
        else
        {
          for (size_t x=min_x; x<=max_x; ++x)
          {
            alpha = (F(x) - b.x) / (a.x - b.x);
            const size_t y = size_t(b.y + alpha * (a.y - b.y) + F(0.5));
            I(x, y) = ImageFloat(F(1) - alpha) * c2 + ImageFloat(alpha) * c1;
          }
        }

      }
      else
      {

        if (a.y < b.y)
        {
          for (size_t y=min_y; y<=max_y; ++y)
          {
            alpha = (F(y) - a.y) / (b.y - a.y);
            const size_t x = size_t(a.x + alpha * (b.x - a.x) + F(0.5));
            I(x, y) = ImageFloat(F(1) - alpha) * c1 + ImageFloat(alpha) * c2;
          }
        }
        else
        {
          for (size_t y=min_y; y<=max_y; ++y)
          {
            alpha = (F(y) - b.y) / (a.y - b.y);
            const size_t x = size_t(b.x + alpha * (a.x - b.x) + F(0.5));
            I(x, y) = ImageFloat(F(1) - alpha) * c2 + ImageFloat(alpha) * c1;
          }
        }

      }

    }

    /**
     * Fill an image of type T with an interpolated triangle. Uses the float type
     * F. Beware: using double actually makes better triangles, see the testcase.
     **/
    template <class F, class Pixel>
    static inline
    void fill_triangle(gravis::tImage<Pixel> &I,
                       const gravis::t2Vector<F> &a,  const gravis::t2Vector<F> &b,  const gravis::t2Vector<F> &c,
                       const Pixel&               c1, const Pixel&               c2, const Pixel&               c3)
    {

      typedef typename tImageTraits<Pixel>::Float_t      ImageFloat;

      // Divisor
      const F det = a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]);

      if (det == F(0)) // The triangle is singular, it has no area
        return;


      // Bounding Box
      const F min_x_f = std::max<F>(F(0),             std::min(std::min(a[0], b[0]), c[0]));
      const F max_x_f = std::min<F>(F(I.cols()-1),    std::max(std::max(a[0], b[0]), c[0]));
      const F min_y_f = std::max<F>(F(0),             std::min(std::min(a[1], b[1]), c[1]));
      const F max_y_f = std::min<F>(F(I.rows()-1),    std::max(std::max(a[1], b[1]), c[1]));

      // Bounding Box
      const int min_x = int( ceil(  min_x_f ) );
      const int max_x = int( floor( max_x_f ) );
      const int min_y = int( ceil(  min_y_f ) );
      const int max_y = int( floor( max_y_f ) );

      t3Vector<F> lambda;
      t3Vector<F> lambdaWorld;
      // Loop over bounding box
      for (int x = min_x; x <= max_x; ++x)
      {
        bool found = false; // Good for larger triangles, but may slow it down for small triangles. This is anyhow a relatively slow method, but obviously correct.
        for (int y = min_y; y <= max_y; ++y)
        {

          // Gets Barycentric Coordinates in Screen space
          lambda[0] = (F(x) * (b[1] - c[1]) + b[0] * (c[1] - F(y)) + c[0] * (F(y) - b[1])) / det;
          lambda[1] = (a[0] * (F(y) - c[1]) + F(x) * (c[1] - a[1]) + c[0] * (a[1] - F(y))) / det;
          lambda[2] = F(1) - lambda[0] - lambda[1];

          // Test if inside triangle
          if ((F(0) <= lambda[0]) && (F(0) <= lambda[1]) && (F(0) <= lambda[2]))
          {
            found = true;

            I(x,y) = ImageFloat(lambda.x) * c1 + ImageFloat(lambda.y) * c2 + ImageFloat(lambda.z) * c3;

          }
          else if (found)
            break;
        }
      }

    }

    template <class F, class Pixel>
    static inline
    void draw_circle(gravis::tImage<Pixel> &I,
                     const gravis::t2Vector<F> &center, const F& radius,
                     const Pixel& c)
    {
      const F pi=3.14159265358979323846264338327950288419716939937510;
      const int w=I.cols();
      const int h=I.rows();
      for (F a=0; a<=0.25*pi; a+=1.0/(2.0*pi*radius))
      {
        const int ds=int( radius*sin(a) );
        const int dc=int( radius*cos(a) );
        const int x1 = int(center.x + ds);
        const int y1 = int(center.y + dc);
        const int x2 = int(center.x - ds);
        const int y2 = int(center.y - dc);
        const int x3 = int(center.x + dc);
        const int y3 = int(center.y + ds);
        const int x4 = int(center.x - dc);
        const int y4 = int(center.y - ds);

        if (0<=x1 && x1<w && 0<=y1 && y1<h) I(x1, y1) = c;
        if (0<=x2 && x2<w && 0<=y2 && y2<h) I(x2, y2) = c;
        if (0<=x1 && x1<w && 0<=y2 && y2<h) I(x1, y2) = c;
        if (0<=x2 && x2<w && 0<=y1 && y1<h) I(x2, y1) = c;

        if (0<=x3 && x3<w && 0<=y3 && y3<h) I(x3, y3) = c;
        if (0<=x4 && x4<w && 0<=y4 && y4<h) I(x4, y4) = c;
        if (0<=x3 && x3<w && 0<=y4 && y4<h) I(x3, y4) = c;
        if (0<=x4 && x4<w && 0<=y3 && y3<h) I(x4, y3) = c;
      }
    }

    /*!
     * Combine two images into one. Inserts img2 into img1 at the position given
     * by col, row.  Parts of img2 that would be inserted outside of img1 are
     * clipped.
     */
    template <class T1, class T2>
    tImage<T1> inset(const tImage<T1> &img1, const tImage<T2> &img2, const int col, const int row)
    {

      int min_x = -std::min(0, col);
      int min_y = -std::min(0, row);
      int max_x =  std::min((int)img2.cols(), (int)img1.cols()-col);
      int max_y =  std::min((int)img2.rows(), (int)img1.rows()-row);

      tImage<T1> lhs = img1.clone();

      for (int y=min_y; y<max_y; y++)
        for (int x=min_x; x<max_x; x++)
          lhs(x+col,y+row) = img2(x,y);

      return lhs;
    }

  }

}

#endif
