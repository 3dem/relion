/******************************************************************************
 **        Title: Raw array in and output
 **  Description: Input/Output for tArrays.
 **
 **       Author: Brian Amberg, 2006
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#ifndef __GRAVIS_IO_ARRAY_RAW__
#define __GRAVIS_IO_ARRAY_RAW__

#include "../../tArray.h"
#include <fstream>
#include <string>

namespace gravis
{

  namespace io
  {

    /**
     * Raw dump of array contents
     **/
    class ArrayRaw
    {

      public:

        /**
         * Load an Array from a raw file.
         * The only checking possible is for consistent filesize
         **/
        template<class T>
        static
        void load(tArray<T> &out, const std::string& filename)
        {
          std::ifstream is(filename.c_str(), std::ios_base::binary);
          // get length of file:
          is.seekg(0, std::ios::end);
          size_t length = is.tellg();
          is.seekg(0, std::ios::beg);

          if (length / sizeof(T) * sizeof(T) != length)
            GRAVIS_THROW3(Exception, "Invalid array file. The length does not fit.", filename);

          length = length / sizeof(T);

          // Load data
          out.setSize(length);
          is.read((char*)(&out[0]), sizeof(T) * length);
        }

        /**
         * Dump an array to a raw file
         **/
        template<class T>
        static
        void save(const std::string& filename, const tArray<T> &in)
        {
          std::ofstream os(filename.c_str(),  std::ios_base::binary);
          os.write((char*)(&in[0]), sizeof(T) * in.size());
        }

    };

  }
}
#endif
