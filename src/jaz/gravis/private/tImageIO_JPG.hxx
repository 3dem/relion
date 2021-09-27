/******************************************************************************
 **        Title: tImageIO_JPG.hxx
 **  Description: Implements reader/writer for JPG image format.
 **
 **       Author: Jean-Sebastien Pierrard, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#ifndef __GRAVIS__TIMAGE_IO_JPG__
#define __GRAVIS__TIMAGE_IO_JPG__
#ifdef JPEG_FOUND

#include <iostream>
#include <fstream>

#include <jpeglib.h>
#include <setjmp.h>
#include <jerror.h>

namespace gravis
{
  namespace priv
  {

    template <class T>
    class JPGImageReader
    {
        struct gravis_jpg_error_mgr
        {
          struct jpeg_error_mgr pub;	/* "public" fields */

          jmp_buf setjmp_buffer;	/* for return to caller */
        };

        /*
         * Here's the routine that will replace the standard error_exit method:
         */
        static
        void gravis_jpg_error_quiet (j_common_ptr cinfo)
        {
          gravis_jpg_error_mgr* myerr = (gravis_jpg_error_mgr*) cinfo->err;
          /* Return control to the setjmp point */
          longjmp(myerr->setjmp_buffer, 1);
        }
        /*
         * Here's the routine that will replace the standard error_exit method:
         */
        static
        void gravis_jpg_error_exit (j_common_ptr cinfo)
        {
          gravis_jpg_error_mgr* myerr = (gravis_jpg_error_mgr*) cinfo->err;

          /* Always display the message. */
          /* We could postpone this until after returning, if we chose. */
          (*cinfo->err->output_message) (cinfo);

          /* Return control to the setjmp point */
          longjmp(myerr->setjmp_buffer, 1);
        }

      public:
        JPGImageReader () {};
        ~JPGImageReader () {};

        void read (tImage<T>&, const std::string&);

        static
        bool canHandle (const std::string& filename);
    };


    template <class T>
    class JPGImageWriter
    {
        struct gravis_jpg_error_mgr
        {
          struct jpeg_error_mgr pub;	/* "public" fields */

          jmp_buf setjmp_buffer;	/* for return to caller */
        };

        /*
         * Here's the routine that will replace the standard error_exit method:
         */
        static
        void gravis_jpg_error_quiet (j_common_ptr cinfo)
        {
          gravis_jpg_error_mgr* myerr = (gravis_jpg_error_mgr*) cinfo->err;
          /* Return control to the setjmp point */
          longjmp(myerr->setjmp_buffer, 1);
        }
        /*
         * Here's the routine that will replace the standard error_exit method:
         */
        static
        void gravis_jpg_error_exit (j_common_ptr cinfo)
        {
          gravis_jpg_error_mgr* myerr = (gravis_jpg_error_mgr*) cinfo->err;

          /* Always display the message. */
          /* We could postpone this until after returning, if we chose. */
          (*cinfo->err->output_message) (cinfo);

          /* Return control to the setjmp point */
          longjmp(myerr->setjmp_buffer, 1);
        }

      public:
        JPGImageWriter () {};
        ~JPGImageWriter () {};

        void write (const tImage<T>& image, const std::string& filename, int quality=100);
    };


    template <class T>
    inline
    bool JPGImageReader<T>::canHandle (const std::string& filename)
    {
      /* This struct contains the JPEG decompression parameters and pointers to
       * working space (which is allocated as needed by the JPEG library).
       */
      struct jpeg_decompress_struct cinfo;

      FILE* infile;		/* source file */
      if ((infile = fopen(filename.c_str(), "rb")) == NULL)
      {
		GRAVIS_THROW3(gravis::Exception, "Could not open file", filename);
      }
      /* We set up the normal JPEG error routines, then override error_exit. */
      struct gravis_jpg_error_mgr jerr;
      cinfo.err = jpeg_std_error(&jerr.pub);
      jerr.pub.error_exit = gravis_jpg_error_quiet;

      /* Establish the setjmp return context for gravis_jpg_error_exit to use. */
      if (setjmp(jerr.setjmp_buffer))
      {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return false;
      }
      jpeg_create_decompress(&cinfo);
      jpeg_stdio_src(&cinfo, infile);
      (void) jpeg_read_header(&cinfo, TRUE);
      return true;
    }


    template <class T>
    inline
    void JPGImageReader<T>::read (tImage<T>& image, const std::string& filename)
    {

      FILE* infile = fopen(filename.c_str(), "rb");
      if (infile == NULL)
      {
        GRAVIS_THROW3(Exception, "Unable to open file: ", filename);
      }

      /* This struct contains the JPEG decompression parameters and pointers to
       * working space (which is allocated as needed by the JPEG library).
       */
      struct jpeg_decompress_struct cinfo;
      /* We use our private extension JPEG error handler.
       * Note that this struct must live as long as the main JPEG parameter
       * struct, to avoid dangling-pointer problems.
       */
      struct gravis_jpg_error_mgr jerr;
      /* More stuff */
      JSAMPARRAY buffer;		/* Output row buffer */
      int row_stride;		/* physical row width in output buffer */

      /* Step 1: allocate and initialize JPEG decompression object */

      /* We set up the normal JPEG error routines, then override error_exit. */
      cinfo.err = jpeg_std_error(&jerr.pub);
      jerr.pub.error_exit = gravis_jpg_error_exit;
      /* Establish the setjmp return context for my_error_exit to use. */
      if (setjmp(jerr.setjmp_buffer))
      {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
		GRAVIS_THROW3(gravis::Exception, "Could not read jpeg file", filename);
      }
      /* Now we can initialize the JPEG decompression object. */
      jpeg_create_decompress(&cinfo);

      /* Step 2: specify data source (eg, a file) */

      jpeg_stdio_src(&cinfo, infile);

      /* Step 3: read file parameters with jpeg_read_header() */

      (void) jpeg_read_header(&cinfo, TRUE);
      /* We can ignore the return value from jpeg_read_header since
       *   (a) suspension is not possible with the stdio data source, and
       *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
       * See libjpeg.doc for more info.
       */

      /* Step 4: set parameters for decompression */

      /* In this example, we don't need to change any of the defaults set by
       * jpeg_read_header(), so we do nothing here.
       */

      /* Step 5: Start decompressor */

      (void) jpeg_start_decompress(&cinfo);
      /* We can ignore the return value since suspension is not possible
       * with the stdio data source.
       */

      /* We may need to do some setup of our own at this point before reading
       * the data.  After jpeg_start_decompress() we have the correct scaled
       * output image dimensions available, as well as the output colormap
       * if we asked for color quantization.
       * In this example, we need to make an output work buffer of the right size.
       */
      image.setSize(cinfo.output_width, cinfo.output_height);
      /* JSAMPLEs per row in output buffer */
      row_stride = cinfo.output_width * cinfo.output_components;
      /* Make a one-row-high sample array that will go away when done with image */
      buffer = (*cinfo.mem->alloc_sarray)
               ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

      /* Step 6: while (scan lines remain to be read) */
      /*           jpeg_read_scanlines(...); */

      /* Here we use the library's state variable cinfo.output_scanline as the
       * loop counter, so that we don't have to keep track ourselves.
       */
      while (cinfo.output_scanline < cinfo.output_height)
      {
        /* jpeg_read_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could ask for
         * more than one scanline at a time if that's more convenient.
         */
        (void) jpeg_read_scanlines(&cinfo, buffer, 1);
        /* Assume put_scanline_someplace wants a pointer and sample count. */
        std::cout << cinfo.output_scanline << std::endl;
        for (size_t x=0; x<cinfo.output_width; ++x)
        {
          const gravis::tRGB<unsigned char> sample(
            buffer[0][3*x], buffer[0][3*x+1], buffer[0][3*x+2]);
          T& pixel = image(x,cinfo.output_scanline-1);
          pixelTypeConverter(sample, pixel);
        }
      }

      /* Step 7: Finish decompression */

      (void) jpeg_finish_decompress(&cinfo);
      /* We can ignore the return value since suspension is not possible
       * with the stdio data source.
       */

      /* Step 8: Release JPEG decompression object */

      /* This is an important step since it will release a good deal of memory. */
      jpeg_destroy_decompress(&cinfo);

      /* After finish_decompress, we can close the input file.
       * Here we postpone it until after no more JPEG errors are possible,
       * so as to simplify the setjmp error logic above.  (Actually, I don't
       * think that jpeg_destroy can do an error exit, but why assume anything...)
       */
      fclose(infile);

      /* At this point you may want to check to see whether any corrupt-data
       * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
       */
      if (jerr.pub.num_warnings != 0)
        GRAVIS_THROW3(gravis::Exception, "Corrupted data in jpeg file", filename);

      image.setName(std::string(filename));

    }




    template <class T>
    inline
    void JPGImageWriter<T>::write (const tImage<T>& image, const std::string& filename, int quality)
    {

      tImageConverter<tImage<T>, tImage< tRGB<unsigned char> > > convert;
      gravis::tImage< tRGB< unsigned char > > rgbimage = convert(image);

      FILE* outfile = fopen(filename.c_str(), "wb");
      if (outfile == 0)
      {
        GRAVIS_THROW3(Exception, "Unable to open/create image file: ", filename);
      }

      if (image.cols() == 0 || image.rows()==0)
        GRAVIS_THROW3(Exception, "Can not write an empty image", filename);

      /* This struct contains the JPEG compression parameters and pointers to
       * working space (which is allocated as needed by the JPEG library).
       * It is possible to have several such structures, representing multiple
       * compression/decompression processes, in existence at once.  We refer
       * to any one struct (and its associated working data) as a "JPEG object".
       */
      struct jpeg_compress_struct cinfo;
      /* This struct represents a JPEG error handler.  It is declared separately
       * because applications often want to supply a specialized error handler
       * (see the second half of this file for an example).  But here we just
       * take the easy way out and use the standard error handler, which will
       * print a message on stderr and call exit() if compression fails.
       * Note that this struct must live as long as the main JPEG parameter
       * struct, to avoid dangling-pointer problems.
       */
      struct gravis_jpg_error_mgr jerr;
      cinfo.err = jpeg_std_error(&jerr.pub);
      jerr.pub.error_exit = gravis_jpg_error_exit;
      /* Establish the setjmp return context for my_error_exit to use. */
      if (setjmp(jerr.setjmp_buffer))
      {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_compress(&cinfo);
        fclose(outfile);
        return;
      }
      /* More stuff */
      JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
      int row_stride;		/* physical row width in image buffer */

      /* Step 1: allocate and initialize JPEG compression object */

      /* We have to set up the error handler first, in case the initialization
       * step fails.  (Unlikely, but it could happen if you are out of memory.)
       * This routine fills in the contents of struct jerr, and returns jerr's
       * address which we place into the link field in cinfo.
       */
      cinfo.err = jpeg_std_error(&jerr.pub);
      /* Now we can initialize the JPEG compression object. */
      jpeg_create_compress(&cinfo);

      /* Step 2: specify data destination (eg, a file) */
      jpeg_stdio_dest(&cinfo, outfile);

      /* Step 3: set parameters for compression */

      /* First we supply a description of the input image.
       * Four fields of the cinfo struct must be filled in:
       */
      cinfo.image_width = image.cols(); 	/* image width and height, in pixels */
      cinfo.image_height = image.rows();
      cinfo.input_components = 3;		/* # of color components per pixel */
      cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
      /* Now use the library's routine to set default compression parameters.
       * (You must set at least cinfo.in_color_space before calling this,
       * since the defaults depend on the source color space.)
       */
      jpeg_set_defaults(&cinfo);
      /* Now you can set any non-default parameters you wish to.
       * Here we just illustrate the use of quality (quantization table) scaling:
       */
      jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

      /* Step 4: Start compressor */

      /* TRUE ensures that we will write a complete interchange-JPEG file.
       * Pass TRUE unless you are very sure of what you're doing.
       */
      jpeg_start_compress(&cinfo, TRUE);

      /* Step 5: while (scan lines remain to be written) */
      /*           jpeg_write_scanlines(...); */

      /* Here we use the library's state variable cinfo.next_scanline as the
       * loop counter, so that we don't have to keep track ourselves.
       * To keep things simple, we pass one scanline per call; you can pass
       * more if you wish, though.
       */
      row_stride = image.cols() * 3;	/* JSAMPLEs per row in image_buffer */

      while (cinfo.next_scanline < cinfo.image_height)
      {
        /* jpeg_write_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could pass
         * more than one scanline at a time if that's more convenient.
         */
        row_pointer[0] = &rgbimage(0,cinfo.next_scanline)[0];
        (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
      }

      /* Step 6: Finish compression */

      jpeg_finish_compress(&cinfo);
      /* After finish_compress, we can close the output file. */
      fclose(outfile);

      /* Step 7: release JPEG compression object */

      /* This is an important step since it will release a good deal of memory. */
      jpeg_destroy_compress(&cinfo);

    }

  } /* Close Namespace "priv" */
} /* Close Namespace "gravis" */
#endif
#endif
