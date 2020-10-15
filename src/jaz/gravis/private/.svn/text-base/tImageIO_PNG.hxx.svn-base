/******************************************************************************
**        Title: tImageIO_PNG.hxx
**  Description: Implements reader/writer for PNG image format.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#include <iostream>
#include <fstream>

#include <png.h>

namespace gravis
{
  namespace priv
  {

    template <class T>
    class PNGImageReader
    {
      public:
        PNGImageReader ();
        ~PNGImageReader ();

        void read (tImage<T>&, const char*);

        static
        bool canHandle (const char*);
    };


    template <class T>
    class PNGImageWriter
    {
      public:
        PNGImageWriter ();
        ~PNGImageWriter ();

        void write (const tImage<T>& image, const char* filename);
    };


    template <class T>
    inline
    PNGImageReader<T>::PNGImageReader () { }


    template <class T>
    inline
    PNGImageReader<T>::~PNGImageReader () { }


    template <class T>
    inline
    bool PNGImageReader<T>::canHandle (const char* header)
    {
      bool is_png = !png_sig_cmp((png_byte*)header, 0, 4);
      return is_png;
    }


    template <class T>
    inline
    void PNGImageReader<T>::read (tImage<T>& image, const char* filename)
    {

      FILE* pngfilep = fopen(filename, "rb");
      if (pngfilep == 0)
      {
        GRAVIS_THROW3(Exception, "Unable to open file: ", filename);
      }

      png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, this, 0, 0);
      if (png_ptr == 0)
      {
        fclose(pngfilep);
        GRAVIS_THROW3(Exception, "PNG Read Error: ", filename);
      }


      png_infop info_ptr = png_create_info_struct(png_ptr);
      if (info_ptr == 0)
      {
        fclose(pngfilep);
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        GRAVIS_THROW3(Exception, "PNG Read Error: ", filename);
      }


      png_init_io(png_ptr, pngfilep);
      png_read_info(png_ptr, info_ptr);


      png_uint_32 width, height;
      int         depth, color_type, interlace_type;

      png_get_IHDR(
        png_ptr, info_ptr,
        &width, &height, &depth, &color_type, &interlace_type,
        NULL, NULL
      );


      // Reduce 16 bit images to 8 bit
      png_set_strip_16(png_ptr);

      // Expand paletted images to RGB
      if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png_ptr);

      // Expand grayscale images to the full 8 bits from 1, 2, or 4 bits/pixel
      if ((color_type == PNG_COLOR_TYPE_GRAY) && (depth < 8))
        png_set_expand_gray_1_2_4_to_8(png_ptr);

      // Expang grayscale image to RGB(A),
      // This should definitely be optimzed !!
      if (color_type == PNG_COLOR_TYPE_GRAY)
      {
        png_set_gray_to_rgb(png_ptr);
      }

      // /* ???
      // Add filler (or alpha) byte (before/after each RGB triplet)
      png_set_filler(png_ptr, 0xff, PNG_FILLER_AFTER);
      // ??? */
      //
      png_set_interlace_handling(png_ptr);

      png_read_update_info(png_ptr, info_ptr); // ???

      // std::cout << "Image: " << width << " x " << height << "(" << depth << "bpp)\n";
      // std::cout << "Bytes/Row: " <<  png_get_rowbytes(png_ptr, info_ptr) << "\n";

      png_bytep* row_pointers = new png_bytep[height];

      for (unsigned int row=0; row<height; row++)
      {
        row_pointers[row] = (png_bytep)png_malloc(png_ptr, png_get_rowbytes(png_ptr, info_ptr));
      }

      // Read entire image
      png_read_image(png_ptr, row_pointers);

      // Read remains
      png_read_end(png_ptr, info_ptr);

      image.setSize(width, height);
      image.setName(std::string(filename));

      for (unsigned int row=0; row<height; ++row)
      {
        T*             tgt_data = image.data(0, row);
        unsigned char* src_data = row_pointers[row];

        for (unsigned int col=0; col<width; ++col)
        {

          tRGBA<unsigned char> pixel(src_data[0], src_data[1], src_data[2], src_data[3]);
          pixelTypeConverter(pixel, *tgt_data);

          tgt_data++;
          src_data += 4;
        }
      }

      // Free PNG-image
      for (unsigned int row=0; row<height; ++row)
      {
        png_free(png_ptr, row_pointers[row]);
      }

      delete[] row_pointers;

      // Cleanup png structures and open file
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(pngfilep);
    }



    template <class T>
    inline
    PNGImageWriter<T>::PNGImageWriter ()
    {
    }

    template <class T>
    inline
    PNGImageWriter<T>::~PNGImageWriter ()
    {
    }


    template <class T>
    inline
    void PNGImageWriter<T>::write (const tImage<T>& image, const char* filename)
    {

      FILE* pngfilep = fopen(filename, "wb");
      if (pngfilep == 0)
      {
        GRAVIS_THROW3(Exception, "Unable to open/create image file: ", filename);
      }

      png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, this, 0, 0);
      if (png_ptr == 0)
      {
        fclose(pngfilep);
        GRAVIS_THROW3(Exception, "PNG(internal) write error: ", filename);
      }

      png_infop info_ptr = png_create_info_struct(png_ptr);
      if (info_ptr == 0)
      {
        fclose(pngfilep);
        png_destroy_write_struct(&png_ptr,  NULL);
        GRAVIS_THROW3(Exception, "PNG(internal) write error: ", filename);
      }

      png_init_io(png_ptr, pngfilep);

      int color_type;

      int nof_cc = tImage_Traits<tImage<T> >::nofComponents();

      switch (nof_cc)
      {
        case 1 :
          color_type = PNG_COLOR_TYPE_GRAY;
          break;
        case 3 :
          color_type = PNG_COLOR_TYPE_RGB;
          break;
        case 4 :
          color_type = PNG_COLOR_TYPE_RGB_ALPHA;
          break;

        default :
        {
          std::cerr << "Unhandled number of color components \n" ;
          exit(1);
        }
      }

      if (image.cols() == 0 || image.rows()==0)
        GRAVIS_THROW3(Exception, "Can not write an empty image", filename);

      png_set_IHDR(
        png_ptr, info_ptr, image.cols(), image.rows(), 8/*depth*/,
        color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE
      );


      png_write_info(png_ptr, info_ptr);

      switch (nof_cc)
      {
        case 1 :
        {
          tImageConverter<tImage<T>, tImage<unsigned char> > convert;
          tImage<unsigned char> out(image.cols(), image.rows());
          out = convert(image);

          png_bytep* row_pointers = new png_bytep[image.rows()];
          for (unsigned int i=0; i<image.rows(); ++i)
          {
            row_pointers[i] = (unsigned char*)out.data(0, i);
          }

          png_write_image(png_ptr, row_pointers);
          delete[] row_pointers;
          break;
        }

        case 3 :
        {
          tImageConverter<tImage<T>, tImage<bRGB> > convert;
          tImage<bRGB> out(image.cols(), image.rows());
          out = convert(image);

          png_bytep* row_pointers = new png_bytep[image.rows()];
          for (unsigned int i=0; i<image.rows(); ++i)
          {
            row_pointers[i] = (unsigned char*)out.data(0, i);
          }

          png_write_image(png_ptr, row_pointers);
          delete[] row_pointers;
          break;
        }

        case 4 :
        {
          tImageConverter<tImage<T>, tImage<bRGBA> > convert;
          tImage<bRGBA> out(image.cols(), image.rows());
          out = convert(image);

          png_bytep* row_pointers = new png_bytep[image.rows()];
          for (unsigned int i=0; i<image.rows(); ++i)
          {
            row_pointers[i] = (unsigned char*)out.data(0, i);
          }

          png_write_image(png_ptr, row_pointers);
          delete[] row_pointers;
          break;
        }

        default :
        {
          std::cerr << "Unhandled number of color components \n" ;
          exit(1);
        }
      }

      png_write_end(png_ptr, info_ptr);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      fclose(pngfilep);
    }

  } /* Close Namespace "priv" */
} /* Close Namespace "gravis" */
