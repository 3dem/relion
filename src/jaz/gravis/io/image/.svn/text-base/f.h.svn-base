/******************************************************************************
**        Title: gravis/io/image/f.h
**  Description: Implements reader/writer for the .f (vector) file format
**
**       Author: Pascal Paysan
**               Brian Amberg
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef __GRAVIS__IO_IMAGE_F__
#define __GRAVIS__IO_IMAGE_F__
#include <iostream>
#include <sstream>
#include <vector>
#include "../../t2Vector.h"
#include "../../t3Vector.h"
#include "../../t4Vector.h"
#include <vector>
#include <zlib.h>

namespace gravis
{
  namespace io
  {

    static const unsigned int FI_VERSION = 1u;
    static const float FI_MAGIC_NUMBER = 3.141592653f*3.141592653f;

    /**
     *
     * saves and load .f displacement fields into images.
     *
     **/

    class ImageF
    {
      private:
        typedef struct FImageHeader
        {
          char name[128];
          unsigned int width;
          unsigned int height;
          unsigned int depth;
          unsigned int channels;
          unsigned int version; //20 byte
          float magicNo; //4 byte
          char reserved[360];
        };

        static inline
        void displayHeader(const FImageHeader& fih)
        {
          std::cout << "Name:     " << fih.name << std::endl;
          std::cout << "Width:    " << fih.width << std::endl;
          std::cout << "Height:   " << fih.height << std::endl;
          std::cout << "Depth:    " << fih.depth << std::endl;
          std::cout << "Version:  " << fih.version << std::endl;
          std::cout << "MagicNo:  " << fih.magicNo << std::endl;
          std::cout << "Channels: " << fih.channels << std::endl;
        };

      public:
        template <class T>
        static
        inline
        void load(::gravis::tImage<T> &image, const std::string& filename)
        {
          typedef typename tImageTraits< T >::Scalar_t T_scalar;
          const size_t channels = sizeof(T) / sizeof(T_scalar);

          gzFile fin = gzopen(filename.c_str(), "rb");
          try
          {
            if (0 == fin)
            {
              GRAVIS_THROW3(Exception, "Unable to open file: ", filename.c_str());
            }


            FImageHeader fih;
            gzread(fin, (char*)(&fih),  sizeof(fih));

            if (!fih.magicNo == FI_MAGIC_NUMBER)
              GRAVIS_THROW3(Exception, "Not a .f file: ", filename);

            if(fih.magicNo  != FI_MAGIC_NUMBER)
              GRAVIS_THROW3(Exception, "File is not an FImage: ", filename);

            if(fih.depth != 1)
              GRAVIS_THROW3(Exception, "Unable to open .f file with depth greater one: ", filename);

            if(fih.version > FI_VERSION)
            {
              std::stringstream s;
              s << "Unable to FImage file with version greater " << FI_VERSION << ". This file has FI_VERSION: "  << fih.version;
              GRAVIS_THROW3(Exception, s.str(), filename);
            }

            if(fih.channels != channels)
            {
              std::stringstream s;
              s << "The image datatype has " << channels << " dimensions, while the .f file has "  << fih.channels << " dimensions." <<
                " Can not load this .f file into the image.";
              GRAVIS_THROW3(Exception, s.str(), filename);
            }

            image.setName(fih.name);
            image.setSize(fih.width, fih.height);

            T_scalar* image_data = reinterpret_cast<T_scalar*>(&image[0]);
            std::vector<float> fbuf(fih.width*fih.height*fih.channels);
            gzread(fin, (char*)(&fbuf[0]), sizeof(fbuf[0])*fbuf.size());
            for(size_t i=0; i<fbuf.size(); ++i)
              image_data[i] = T_scalar(fbuf[i]);

          }
          catch (const char* e)
          {
            gzclose(fin);
            throw(e);
          }
          catch (gravis::Exception& e)
          {
            gzclose(fin);
            throw(e);
          }
          gzclose(fin);
        }

        template <class T>
        static
        inline
        void save(const std::string& filename, const ::gravis::tImage<T> &image)
        {
          typedef typename tImageTraits< T >::Scalar_t T_scalar;
          const size_t channels = sizeof(T) / sizeof(T_scalar);

          std::ofstream of(filename.c_str(), std::ofstream::binary);
          if (!of.good())
          {
            GRAVIS_THROW3(Exception, "Unable to open/create .f file: ", filename);
          }

          FImageHeader fih;
          strncpy(fih.name, image.name().c_str(), 128);
          fih.width    = image.cols();
          fih.height   = image.rows();
          fih.depth    = 1;
          fih.version  = FI_VERSION;
          fih.magicNo  = FI_MAGIC_NUMBER;
          fih.channels = channels;
          of.write((char*)(&fih),sizeof(FImageHeader));

          const T_scalar* image_data = reinterpret_cast<const T_scalar*>(&image[0]);
          std::vector<float> fbuf(fih.width*fih.height*fih.channels);
          for(size_t i=0; i<fbuf.size(); ++i)
            fbuf[i] = float(image_data[i]);

          of.write((char*)(&fbuf[0]),sizeof(fbuf[0])*fbuf.size());
        }
    };
  }
}
#endif
