/******************************************************************************
**        Title: tImageIO_PNM.hxx
**  Description: Implements reader/writer for PNM image format.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#include <iostream>
#include <fstream>

namespace gravis
{
  namespace priv
  {

    template <class T>
    class PNMImageReader
    {
      public:

        PNMImageReader ();
        ~PNMImageReader ();

        void read (tImage<T>& image, const char* filename);

        static
        bool canHandle (const char* header);

      protected:
        bool garbage;

        bool isWS   (int) const;
        int  getInt (std::istream&);
    };


    template <class T>
    class PNMImageWriter
    {
      public:
        PNMImageWriter ();
        ~PNMImageWriter ();

        void write (const tImage<T>& image, const char* filename);

      public:
        bool raw, fullcolor;
    };



    template <class T>
    inline
    PNMImageReader<T>::PNMImageReader () { }


    template <class T>
    inline
    PNMImageReader<T>::~PNMImageReader () { }


    template <class T>
    inline
    bool PNMImageReader<T>::isWS (int c) const
    {
      return ((c == ' ') || (c == '\r') || (c == '\t') || (c == '\n'));
    }

    template <class T>
    int PNMImageReader<T>::getInt (std::istream& is)
    {
      int c, i = 0;

      c = is.get();
      while (1)
      {

        // Skip comments (lines starting with #)
        if (c == '#')
        {
          while (1)
          {
            c = is.get();
            if ((c == '\n') || (c == EOF)) break;
          }
        }

        if (c == EOF) return EOF;
        if ((c >= '0') && (c <= '9')) break;
        if (!isWS(c)) garbage = true;

        c = is.get();
      }

      while (1)
      {
        i = (i * 10) + (c - '0');
        c = is.get();
        if (c == EOF) return i;
        if ((c < '0') || (c > '9')) break;
      }
      return i;
    }


    // Identify PNM images by file content (header)
    // Search for "P3\n"(ascii) or "P6\n"(binary) RGB
    //            "P2\n"(ascii) or "P5\n"(binary) GRAYSCALE

    template <class T>
    inline
    bool PNMImageReader<T>::canHandle (const char* header)
    {
      if ( (header[0] == 'P') &&
           ((header[1] == '3') || (header[1] == '6') ||
            (header[1] == '2') || (header[1] == '5')) &&
           (header[2] == (char)0x0A) ) return true;
      else return false;
    }


    template <class T>
    inline
    void PNMImageReader<T>::read (tImage<T>& image, const char* filename)
    {

      int  wd, ht, colors, c, counter = 0;
      bool grayscale = true, ascii = true;

      std::ifstream is(filename, std::ios::in | std::ios::binary);
      if (!is.good())
      {
        GRAVIS_THROW3(Exception, "Unable to open file", filename);
      }

      garbage = false;

      c = is.get();    // read 'P'
      c = is.get();    // read format identifier
      if (c == '2')
      {
        grayscale = true;
        ascii = true;
      }
      if (c == '3')
      {
        grayscale = false;
        ascii = true;
      }
      if (c == '5')
      {
        grayscale = true;
        ascii = false;
      }
      if (c == '6')
      {
        grayscale = false;
        ascii = false;
      }

      wd = getInt(is);
      ht = getInt(is);
      colors = getInt(is);

      if (garbage)
      {
        GRAVIS_THROW3(Exception, "Corrupt image header", filename);
      }

      if ((wd <= 0) || (ht <= 0))
      {
        GRAVIS_THROW3(Exception, "Illegal image dimensions", filename);
      }

      if (colors > 255)
      {
        GRAVIS_THROW3(Exception, "Unsupported colorvalues", filename);
      }

      image.setSize(wd, ht);
      image.setName(std::string(filename));

      if (grayscale)
      {
        T* data = image.begin();

        if (ascii)
        {
          for (; data!=image.end(); ++data)
          {
            c = getInt(is);

            if (c >= 0)
            {
              unsigned char value = (unsigned char)c;

              pixelTypeConverter(value, *data);

              ++counter;
            }
            else break;
          }
        }
        else     // not ascii
        {
          unsigned char* bbuf = new unsigned char[image.size()];
          unsigned char* bptr = bbuf;

          is.read((char*)bbuf, image.size());
          counter += is.gcount();

          for (; data!=image.end(); ++data, ++bptr)
          {
            pixelTypeConverter(*bptr, *data);
          }
          delete[] bbuf;
        }

      }
      else     // not grayscale
      {
        T* data = image.begin();

        if (ascii)
        {
          for (; data!=image.end(); ++data)
          {
            bRGB pnm_pixel;

            c = getInt(is);
            if (c >= 0) pnm_pixel.r = (unsigned char)c;
            else break;

            c = getInt(is);
            if (c >= 0) pnm_pixel.g = (unsigned char)c;
            else break;

            c = getInt(is);
            if (c >= 0) pnm_pixel.b = (unsigned char)c;
            else break;

            pixelTypeConverter(pnm_pixel, *data);

            ++counter;
          }
        }
        else     // not grayscale, not ascii
        {

          unsigned char* bbuf = new unsigned char[3*image.cols()];

          for (unsigned int y=0; y<image.rows(); y++)
          {
            unsigned char* bptr = bbuf;

            is.read((char*)bbuf, 3*image.cols());
            counter += is.gcount()/3;

            for (unsigned int i=0; i<image.cols(); ++i)
            {
              bRGB pnm_pixel(bptr[0], bptr[1], bptr[2]);
              pixelTypeConverter(pnm_pixel, *data);

              ++data;
              bptr += 3;
            }
          }

          delete [] bbuf;
        }
      }

      if (garbage)
      {
        std::cerr << "Corrupt image data: " << filename << "\n";
      }
      if (counter != ((int)image.size()))
      {
        std::cerr << "Image appears to be truncated: " << filename << "\n";
      }
    }


    /******************************************************************************
    ******************************************************************************/

    template <class T>
    inline
    PNMImageWriter<T>::PNMImageWriter ()
    {
      raw       = true;
      fullcolor = true;
    }


    template <class T>
    inline
    PNMImageWriter<T>::~PNMImageWriter () { }


    template <class T>
    inline
    void PNMImageWriter<T>::write (const tImage<T>& image, const char* filename)
    {

      std::ofstream os(filename, std::ios::out | std::ios::binary);
      if (!os.good())
      {
        GRAVIS_THROW3(Exception, "Unable to open/create file", filename);
      }

      T* data = image.begin();

      if (fullcolor == false)   // Write GRAY values
      {
        if (raw == false)   // Write 'ascii' format
        {

          os << "P2\n" << image.cols() << " " << image.rows() << "\n255\n";

          int cols_per_line = 0;
          for (; data!=image.end(); ++data)
          {
            if (cols_per_line > 8)
            {
              cols_per_line = 0;
              os << "\n";
            }

            unsigned char value = 0;
            pixelTypeConverter(*data, value);

            os << int(value) << " ";
            ++cols_per_line;
          }

        }
        else     // Write 'raw' format
        {
          os << "P5\n" << image.cols() << " " << image.rows() << "\n255\n";

          unsigned char* bbuf = new unsigned char[image.cols()];

          for (int y=0; y<image.rows(); y++)
          {
            unsigned char* bptr = bbuf;

            // Convert one h-line of pixels into tmp-buffer
            for (int i=0; i<image.cols(); ++i)
            {

              pixelTypeConverter(*data, *bptr);
              ++bptr;
              ++data;
            }

            // Write one h-line of pixels
            os.write((char*)bbuf, image.cols());
          }
          delete[] bbuf;

        }
      }
      else     // Write RGB values
      {
        if (raw == false)   // Write 'ascii' format
        {

          os << "P3\n" << image.cols() << " " << image.rows() << "\n255\n";

          int cols_per_line = 0;
          for (; data!=image.end(); ++data)
          {
            if (cols_per_line > 8)
            {
              cols_per_line = 0;
              os << "\n";
            }

            tRGB<unsigned char> value;
            pixelTypeConverter(*data, value);

            os << int(value.r) << " " << int(value.g) << " " << int(value.b) << " ";
            ++cols_per_line;
          }

        }
        else     // Write 'raw' format
        {

          os << "P6\n" << image.cols() << " " << image.rows() << "\n255\n";

          tRGB<unsigned char>* bbuf = new tRGB<unsigned char>[image.cols()];

          for (int y=0; y<image.rows(); y++)
          {
            tRGB<unsigned char>* bptr = bbuf;

            // Convert one h-line of pixels into tmp-buffer
            for (int i=0; i<image.cols(); ++i)
            {
              pixelTypeConverter(*data, *bptr);
              ++bptr;
              ++data;
            }

            // Write one h-line of pixels
            os.write((char*)bbuf, 3*image.cols());
          }
          delete[] bbuf;

        }
      }

      os.close();
    }


  } /* Close Namespace "priv" */
} /* Close Namespace "gravis" */
