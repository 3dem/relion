/******************************************************************************
 **        Title: Checked access for tImage
 **  Description: Checked Image access
 **
 **       Author: Jean Sebastian Pierrard, 2005
 **               Brian Amberg, 2007
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/

#ifndef __GRAVIS__IMAGE_ACCESS__
#define __GRAVIS__IMAGE_ACCESS__

namespace gravis
{
  namespace image_access
  {

    //! Functor for access behind the image borders
    struct Zero
    {
      template <class T>
      static
      T getPixel (const tImage<T>& image, int x, int y)
      {
        if (x < 0) return T(0);
        else if (x >= (int)image.cols()) return T(0);
        if (y < 0) return T(0);
        else if (y >= (int)image.rows()) return T(0);
        return image(x, y);
      }
    };

    //! Functor for access behind the image borders
    struct Repeat
    {
      template <class T>
      static
      const T& getPixel (const tImage<T>& image, int x, int y)
      {
        if (x < 0) x = 0;
        else if (x >= (int)image.cols()) x = image.cols()-1;
        if (y < 0) y = 0;
        else if (y >= (int)image.rows()) y = image.rows()-1;
        return image(x, y);
      }
    };


    //! Functor for access behind the image borders
    struct Mirror
    {
      template <class T>
      static
      const T& getPixel (const tImage<T>& image, int x, int y)
      {
        if (x < 0) x = 0-x;
        else if (x >= (int)image.cols()) x = 2*(image.cols()-1) - x;
        if (y < 0) y = 0-y;
        else if (y >= (int)image.rows()) y = 2*(image.rows()-1) - y;
        return image(x, y);
      }
    };


    //! Functor for access behind the image borders
    struct Wrap
    {
      template <class T>
      static
      const T& getPixel (const tImage<T>& image, int x, int y)
      {
        if (x < 0) x = image.cols()-1+x;
        else if (x >= (int)image.cols()) x = x-image.cols();
        if (y < 0) y = image.rows()-1+y;
        else if (y >= (int)image.rows()) y = y-image.rows();
        return image(x, y);
      }
    };

  }
}


#endif
