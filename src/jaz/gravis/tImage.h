#ifndef __LIBGRAVIS_T_IMAGE_H__
#define __LIBGRAVIS_T_IMAGE_H__
/******************************************************************************
**        Title: tImage.h
**  Description: Implements two dimensional array with row-major memory layout.
**
**       Author: Jean-Sebastien Pierrard, 2009
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <string>
#include <cstring>

#include "tRGB.h"
#include "tBGR.h"
#include "tRGBA.h"
#include "tRGB_A.h"
#include "tGray_A.h"
#include "tArray.h"
#include <ctype.h>
#include "tImage/traits.h"


/*!
** \file tImage.h
*/

namespace gravis
{
  template <class T>
  class tImage;
}

#include "tImage/access.hxx"
#include "tImage/interpolation.hxx"


namespace gravis
{

  /*!
   ** \class tImage
   ** \brief Implements two dimensional array with row-major memory layout.
   **
   ** This class represents an image of arbitrary pixel type.  TODO
   **
   ** For operations on images look at tImage/operators.h, tImage/???.h.
   */
  template <class T>
  class tImage
  {

      inline static
      bool has_ending(const std::string& filename, const std::string& ending)
      {
        if (filename.size() < ending.size()) return false;
        for (size_t i=0; i<ending.size(); ++i)
          if (tolower(filename[filename.size()-ending.size()+i]) !=
              tolower(ending[i])) return false;
        return true;
      }

    public:


      typedef T   Pixel_t;
      typedef typename tImageTraits<T>::Scalar_t scalar_type;
      typedef T*  iterator;

      tImage ();
      tImage (size_t, size_t, std::string="");
      tImage (size_t, size_t, T const& value );
      tImage (const tImage<T>&);
      tImage& operator=(const tImage<T>&);
      ~tImage ();

      tImage<T>  clone () const;
      tImage<T>& setSize (size_t, size_t);
      tImage<T>& resize (size_t, size_t);
      tImage<T>& setName (std::string);
      tImage<T>& fill (T);

      std::string name () const;
      size_t      cols () const;
      size_t      rows () const;
      size_t      size () const;
      /** Returns the number of components per pixel **/
      size_t      components () const;


      const T& operator() (size_t, size_t) const;
      T& operator() (size_t, size_t);

      /** Returns the component specified by column, row and, component number (channel) **/
      const scalar_type& operator() (size_t, size_t, size_t) const;
      scalar_type& operator() (size_t, size_t, size_t);

      const T& operator [] (size_t) const;
      T& operator [] (size_t);

      /** Returns the component specified by index (as [] operator) and, component number (channel) **/
      const scalar_type& comp(size_t, size_t) const;
      scalar_type& comp(size_t, size_t);

      iterator begin () const;
      iterator end   () const;

      const T* data () const;
      T* data ();

      const T* data (size_t, size_t) const;
      T* data (size_t, size_t);

      void read (const std::string&);

      /**
       * Detect the filetype from the ending.
       **/
      void write(const std::string&) const;
      void writePNM (const std::string&) const;
      void writePNG (const std::string&) const;
#ifdef JPEG_FOUND
      void writeJPG (const std::string&, int quality=100) const;
#endif

      /**
       * Interpolated access to the image
       *
       * Usage
       *
       *   image.interpolate<interpolation::NearestNeighbour>(x, y)
       *   image.interpolate<interpolation::Linear>(x, y)
       *   image.interpolate<interpolation::Cubic>(x, y)
       *
       * See interpolation:: namespace for other methods.
       *
       * Beware:
       *   if using this inside of a templated function or class, you have to write
       *   image.template interpolate<interpolation::Linear>(x, y), which is quite
       *   awfull.
       **/
      template <class InterpolationMethod, class Float>
      inline
      T interpolate(const Float& x, const Float& y) const
      {
        return InterpolationMethod::getPixel(*this, x, y);
      }

      /**
       * Default interpolation mode is Cubic
       **/
      template <class Float>
      inline
      T interpolate(const Float& x, const Float& y) const
      {
        return interpolation::Cubic::getPixel(*this, x, y);
      }


      /**
       * Checked access to the image, with configurable behaviour.
       *
       * Usage
       *
       *   image.access<access::Zero>(x, y)
       *   image.access<access::Repeat>(x, y)
       *   image.access<access::Mirrored>(x, y)
       *   image.access<access::Wrapped>(x, y)
       *
       * Beware:
       *   if using this inside of a templated function or class, you have to write
       *   image.template access<access::Zero>(x, y), which is quite
       *   awfull.
       *
       **/
      template <class AccessMethod>
      inline
      T access(const int& x, const int& y) const
      {
        return AccessMethod::getPixel(*this, x, y);
      }

      /**
       * Default access mode is access::Repeat
       **/
      inline
      T access(const int& x, const int& y) const
      {
        return image_access::Repeat::getPixel(*this, x, y);
      }

      /**
       * tImage Convolution using the access specified access method.
       *
       * The access methods include:
       *  AccessZero
       *  AccessRepeat
       *  AccessWrapped
       *  AccessMirrored
       *
       * Usage:
       *
       *   tImage<fRGBA> result = image.convolve< access::AccessMirrored >(kernel);
       *
       *
       * Beware:
       *   if using this inside of a templated function or class, you have to write
       *   image.template convolve<access::Zero>(kernel), which is quite
       *   awfull.
       **/
      template <class AccessMethod>
      tImage convolve(const tImage< typename tImageTraits<T>::Float_t >& kernel) const
      {
        int klmargin, ktmargin;

        if ((kernel.cols() % 2) == 0)
        {
          klmargin = (kernel.cols() >> 1) - 1;
        }
        else
        {
          klmargin = (kernel.cols() >> 1);
        }

        if ((kernel.rows() % 2) == 0)
        {
          ktmargin = (kernel.rows() >> 1) - 1;
        }
        else
        {
          ktmargin = (kernel.rows() >> 1);
        }

        tImage<T> lhs(cols(), rows());

        for (int r=0; r<(int)rows(); ++r)
        {
          for (int c=0; c<(int)cols(); ++c)
          {
            T sum = T(0);

            for (int ky=0; ky<(int)kernel.rows(); ++ky)
            {
              for (int kx=0; kx<(int)kernel.cols(); ++kx)
              {
                sum += kernel(kx, ky) * access<AccessMethod>(kx-klmargin+c, ky-ktmargin+r);
              }
            }

            lhs(c, r) = sum;
          }
        }

        return lhs;
      }

      /**
       * Default access method is Repeat
       **/
      tImage convolve(const tImage< typename tImageTraits<T>::Float_t >& kernel) const
      {
        return (*this).template convolve<image_access::Repeat>(kernel);
      }

      /** Clamp an image by calling the clamp() method on each element **/
      void clamp()
      {
        for (size_t i=0; i<size(); ++i) (*this)[i].clamp();
      }

    protected:

      std::string p_name;
      size_t wd;
      size_t ht;

      tArray<T>  image;
      tArray<T*> accel;

      iterator    p_begin;
      iterator    p_end;
  };

} /* Close namespace "gravis" */


/******************************************************************************
** tImage<T> implementation
******************************************************************************/

#include "Exception.h"
#include "private/tImageIO.hxx"
#include "private/tImageConverter.hxx"
#include "private/tImageIO_PNM.hxx"
#include "private/tImageIO_PNG.hxx"
#include "private/tImageIO_JPG.hxx"

namespace gravis
{


  /*!
  ** \brief Default constructor.
  */
  template <class T>
  inline
  tImage<T>::tImage ()
    : p_name(""), wd(0), ht(0),
      image(),
      accel(),
      p_begin(),
      p_end()
  {
  }


  /*!
  ** \brief Constructor.
  **
  ** \param width  Set number of columns.
  ** \param height Set number of rows.
  ** \param name   Sets a name for the image (\em optional).
  */
  template <class T>
  inline
  tImage<T>::tImage (size_t width, size_t height, std::string name)
    : p_name(name), wd(width), ht(height),

      // Allocate space for channel data and indexing accelerators
      image(width* height),
      accel(height),
      p_begin(image.data()),
      p_end(image.data()+image.size())
  {

    // Compute pointers to beginning of each line
    for (size_t y=0; y<height; ++y) accel[y] = image.data() + y*width;
  }

  /*!
  ** \brief Constructor.
  **
  ** \param width  Set number of columns.
  ** \param height Set number of rows.
  ** \param value  initialize every pixel to this value.
  */
  template <class T>
  inline
  tImage<T>::tImage (size_t width, size_t height, T const& value)
    : p_name(""), wd(width), ht(height),

      // Allocate space for channel data and indexing accelerators
      image(width* height),
      accel(height),
      p_begin(image.data()),
      p_end(image.data()+image.size())
  {
    // Compute pointers to beginning of each line
    for (size_t y=0; y<height; ++y) accel[y] = image.data() + y*width;
    fill( value );
  }

  /*!
  ** \brief Copy-constructor.
  **
  ** The copy-constructor has reference-semantic, i.e. the image data is not actually
  ** copied. Instead a new handle to the same data is created.
  **
  ** \param rhs
  */
  template <class T>
  inline
  tImage<T>::tImage (const tImage<T>& rhs) :
    p_name (rhs.p_name),
    wd     (rhs.wd),
    ht     (rhs.ht),

    image  (rhs.image),
    accel  (rhs.accel),

    p_begin(rhs.p_begin),
    p_end  (rhs.p_end)
  {
  }

  /*!
  ** \brief Reference Semantic Assignemnt
  **
  ** The assignmment has reference-semantic, i.e. the image data is not actually
  ** copied. Instead a new handle to the same data is created.
  **
  ** \param rhs
  */
  template <class T>
  inline
  tImage<T> &tImage<T>::operator =(const tImage<T>& rhs)
  {
    p_name  = rhs.p_name;
    wd      = rhs.wd;
    ht      = rhs.ht;

    image   = rhs.image;
    accel   = rhs.accel;

    p_begin = rhs.p_begin;
    p_end   = rhs.p_end;
    return *this;
  }


  /*!
  ** \brief Destructor.
  **
  ** Destroy the object(handle). The image data is \em only deleted if no other
  ** instance of this class holds a reference to it.
  */
  template <class T>
  inline
  tImage<T>::~tImage ()
  {
  }


  /*!
  ** \brief Create a deep-copy of the image data.
  **
  ** \return A new tImage<T> object.
  **
  ** \warning This method creates a byte-wise copy of the image data. When
  ** applied to compound types (e.g. T=std::vector<int>) it is very likely to
  ** cause serious problems.
  */
  template <class T>
  tImage<T> tImage<T>::clone () const
  {
    // Allocate new image with same name and dimensions
    tImage result(wd, ht, p_name);

    // Copy the data
    memcpy(result.data(), data(), wd*ht*sizeof(T));

    return result;
  }


  /*!
  ** \brief Resize image.
  **
  ** \param nwd Number of columns in resized image.
  ** \param nht Number of rows in resized image.
  **
  ** \return
  ** \warning The original data is not copied TODO??
  */
  template <class T>
  inline
  tImage<T>& tImage<T>::resize (size_t nwd, size_t nht)
  {
    if ((nwd != wd) || (nht != ht))
      *this = tImage<T>(nwd, nht, p_name);
    return *this;
  }

  /*!
  ** \brief Resize image.
  **
  ** \param nwd Number of columns in resized image.
  ** \param nht Number of rows in resized image.
  **
  ** \return
  ** \warning The original data is not copied TODO??
  */
  template <class T>
  inline
  tImage<T>& tImage<T>::setSize (size_t nwd, size_t nht)
  {
    return resize(nwd, nht);
  }


  template <class T>
  inline
  size_t tImage<T>::rows () const
  {
    return ht;
  }


  template <class T>
  inline
  size_t tImage<T>::cols () const
  {
    return wd;
  }


  template <class T>
  inline
  size_t tImage<T>::size () const
  {
    return image.size();
  }


  template <class T>
  inline
  std::string tImage<T>::name () const
  {
    return p_name;
  }

  template <class T>
  inline
  size_t tImage<T>::components() const
  {
    return tImageTraits<T>::components();
  }

  template <class T>
  inline
  tImage<T>& tImage<T>::setName (std::string name)
  {
    p_name = name;
    return *this;
  }


  template <class T>
  inline
  tImage<T>& tImage<T>::fill (T value)
  {
    image.fill(value);
    return *this;
  }


  template <class T>
  inline
  const T& tImage<T>::operator() (size_t x, size_t y) const
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert( x >= 0 && x < cols() );
    assert( y >= 0 && y < rows() );
#endif
    return (accel[y])[x];
  }


  template <class T>
  inline
  T& tImage<T>::operator() (size_t x, size_t y)
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert( x >= 0 && x < cols() );
    assert( y >= 0 && y < rows() );
#endif
    return (accel[y])[x];
  }

  template <class T>
  inline
  const typename tImage<T>::scalar_type& tImage<T>::operator() (size_t x, size_t y, size_t c) const
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert( x >= 0 && x < cols() );
    assert( y >= 0 && y < rows() );
    assert( c >= 0 && c < tImageTraits<T>::components() );
#endif
    const scalar_type* p = reinterpret_cast<const scalar_type*>(p_begin);
    return p[(y*cols() + x) * tImageTraits<T>::components() + c];
  }


  template <class T>
  inline
  typename tImage<T>::scalar_type& tImage<T>::operator() (size_t x, size_t y, size_t c)
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert( x >= 0 && x < cols() );
    assert( y >= 0 && y < rows() );
    assert( c >= 0 && c < tImageTraits<T>::components() );
#endif
    scalar_type* p = reinterpret_cast<scalar_type*>(p_begin);
    return p[(y*cols() + x) * tImageTraits<T>::components() + c];
  }


  template <class T>
  inline
  const T& tImage<T>::operator[] (size_t n) const
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert((n >= 0) && (n < image.size()));
#endif
    return *(p_begin + n);
  }


  template <class T>
  inline
  T& tImage<T>::operator[] (size_t n)
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert((n >= 0) && (n < image.size()));
#endif
    return *(p_begin + n);
  }

  template <class T>
  inline
  const typename tImage<T>::scalar_type& tImage<T>::comp(size_t n, size_t c) const
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert((n >= 0) && (n < image.size()));
    assert( c >= 0 && c < tImageTraits<T>::components() );
#endif
    const scalar_type* p = reinterpret_cast<const scalar_type*>(p_begin);
    return *(p + n*tImageTraits<T>::components() + c);
  }


  template <class T>
  inline
  typename tImage<T>::scalar_type& tImage<T>::comp(size_t n, size_t c)
  {
#ifdef _GRAVIS_DEBUG_RANGECHECKING_
    assert((n >= 0) && (n < image.size()));
    assert( c >= 0 && c < tImageTraits<T>::components() );
#endif
    scalar_type* p = reinterpret_cast<scalar_type*>(p_begin);
    return *(p + n*tImageTraits<T>::components() + c);
  }


  template <class T>
  inline
  typename tImage<T>::iterator tImage<T>::begin () const
  {
    return p_begin;
  }


  template <class T>
  inline
  typename tImage<T>::iterator tImage<T>::end () const
  {
    return p_end;
  }


  template <class T>
  inline
  const T* tImage<T>::data () const
  {
    return image.data();
  }


  template <class T>
  inline
  T* tImage<T>::data ()
  {
    return image.data();
  }


  template <class T>
  inline
  const T* tImage<T>::data (size_t x, size_t y) const
  {
    return accel[y] + x;
  }


  template <class T>
  inline
  T* tImage<T>::data (size_t x, size_t y)
  {
    return accel[y] + x;
  }


  template <class T>
  inline
  void tImage<T>::read (const std::string& filename)
  {
#ifdef JPEG_FOUND
    if (priv::JPGImageReader<T>::canHandle(filename))
    {
      priv::JPGImageReader<T> reader;
      reader.read(*this,  filename);
      return;
    }
#endif

    char header[512];

    std::ifstream is(filename.c_str(), std::ios::in | std::ios::binary);
    if (!is.good())
    {
      GRAVIS_THROW3(Exception, "Unable to open file", filename);
    }

    is.read(&header[0], sizeof(header));
    is.close();

    if (priv::PNMImageReader<T>::canHandle(header))
    {
      priv::PNMImageReader<T> reader;
      reader.read(*this, filename.c_str());
      return;
    }

    if (priv::PNGImageReader<T>::canHandle(header))
    {
      priv::PNGImageReader<T> reader;
      reader.read(*this,  filename.c_str());
      return;
    }

    GRAVIS_THROW3(gravis::Exception, "Can't handle this file.", filename);
  }

  template <class T>
  inline
  void tImage<T>::write(const std::string& filename) const
  {
    if (has_ending(filename, "jpg") || has_ending(filename, "jpeg"))
    {
#ifdef JPEG_FOUND
     writeJPG(filename);
#else
     GRAVIS_THROW3(gravis::Exception, "libjpeg was not linked during compilation so we cannot read: ", filename);
#endif
    }
    else if (has_ending(filename, "png"))
      writePNG(filename);
    /*else if (has_ending(filename, "pnm"))
      writePNM(filename);*/
    else
      GRAVIS_THROW3(gravis::Exception, "Could not determine filetype from filename: ", filename);
  }

  template <class T>
  inline
  void tImage<T>::writePNM (const std::string& filename) const
  {
    priv::PNMImageWriter<T> writer;
    writer.write(*this, filename.c_str());
  }

  template <class T>
  inline
  void tImage<T>::writePNG (const std::string& filename) const
  {
    priv::PNGImageWriter<T> writer;
    writer.write(*this, filename.c_str());
  }

#ifdef JPEG_FOUND
  template <class T>
  inline
  void tImage<T>::writeJPG (const std::string& filename, int quality) const
  {
    priv::JPGImageWriter<T> writer;
    writer.write(*this, filename.c_str(), quality);
  }
#endif

} /* Close namespace "gravis" */

#endif
