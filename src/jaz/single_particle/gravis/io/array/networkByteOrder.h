/******************************************************************************
 **        Title: tArrayIO.h
 **  Description: Input/Output for tArrays.
 **
 **       Author: Reinhard Knothe / Michael Keller, 2006
 **               Brian Amberg, 2007
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#ifndef _T_ARRAY_IO_NETWORK_BYTE_ORDER_H_
#define _T_ARRAY_IO_NETWORK_BYTE_ORDER_H_

#include <src/jaz/gravis/tArray.h>
#include <src/jaz/gravis/Exception.h>
#include <fstream>
#include <string>

#include <netinet/in.h>

namespace gravis
{

  namespace io
  {

    class ArrayNetworkByteOrder
    {

      public:

        template<class T>
        static
        tArray<T> load(const std::string& filename, size_t count=0)
        {

          tArray<T> ret;
          load(ret,filename,count);
          return ret;
        }

        template<class T>
        static
        void load(tArray<T> &out, const std::string& filename, size_t count=0)
        {
          FILE* in = fopen(filename.c_str(), "rb");
          if (!in)
          {
            GRAVIS_THROW3(Exception, "Unable to open file", filename);
          }
          // determine size
          if (!count)
          {
            fseek(in, 0, SEEK_END);
            count = (int)(ftell(in)/sizeof(T));
            fseek(in, 0, SEEK_SET);
          }
          // write
          int size = sizeof(T)*count/sizeof(int);
          out.resize(count);
          int* outh = (int*)(out.data());
          int* outn = new int[size];
          if (fread(outn, sizeof(int), size, in) != size_t(size))
          {
            fclose(in);
            out.setSize(0);
            delete outn;
            GRAVIS_THROW3(Exception, "Error reading file", filename);
          }
          fclose(in);
          // copy
          for (int i=0; i<size; ++i )
          {
            outh[i] = ntohl(outn[i]);
          }

          delete outn;
        }

        template <class T>
        static
        void save(const std::string& filename, const tArray<T> &in, size_t count=0)
        {
          if (!count) count = in.size();
          FILE* out = fopen(filename.c_str(), "wb");
          if (!out)
          {
            GRAVIS_THROW3(Exception, "Unable to open file", filename);
          }

          int* inh = (int*)(in.data());

          int size = sizeof(T)*count/sizeof(int);
          int* inn = new int[size];
          for (int i=0; i<size; ++i )
          {
            inn[i] = htonl(inh[i]); //
          }

          if (fwrite(inn, sizeof(int), size, out) != size_t(size))
          {
            fclose(out);
            unlink(filename.c_str());
            GRAVIS_THROW3(Exception, "Error writing file", filename);
          }
          fclose(out);
          delete inn;
        }
    };
  }

} // end namespace

#endif
