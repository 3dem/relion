/******************************************************************************
 **        Title: tArrayIO.h
 **  Description: Input/Output for tArrays.
 **
 **       Author: Michael Keller, 2006
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#ifndef _T_LE_H_
#define _T_LE_H_

#include <stdexcept>

#include <gravis/tArray.h>
#include <fstream>
#include <string>

namespace gravis
{
  namespace io
  {

    class ArrayLittleEndian
    {

      public:

        template<class T>
        static tArray<T> load(std::string filename, size_t count=0);
        template<class T>
        static void load(tArray<T>& out, std::string filename, size_t count=0);
        template<class T>
        static void save(std::string filename, const tArray<T> in, size_t count=0);

    }; //  end class

    template<class T>
    tArray<T> ArrayLittleEndian::load(std::string filename, size_t count)
    {

      tArray<T> out;
      load(out,filename,count);
      return out;
    }
    template<class T>
    void ArrayLittleEndian::load(tArray<T>& out, std::string filename, size_t count)
    {
      FILE* in = fopen(filename.c_str(), "rb");
      if (!in)
      {
        throw std::runtime_error("Cannot open '" + filename + "'!");
      }
      // determine size
      if (!count)
      {
        fseek(in, 0, SEEK_END);
        count = (int)(ftell(in)/sizeof(T));
        fseek(in, 0, SEEK_SET);
      }
      // write
      out.setSize(count);
      if (fread(out.data(), sizeof(T), count, in) != count)
      {
        fclose(in);
        out.setSize(0);
        throw std::runtime_error("Error reading '" + filename + "'!");
      }
      fclose(in);
    }

    template<class T>
    void ArrayLittleEndian::save(std::string filename, const tArray<T> in, size_t count)
    {
      if (!count) count = in.size();
      FILE* out = fopen(filename.c_str(), "wb");
      if (!out)
      {
        throw std::runtime_error("Cannot open '" + filename + "'!");
      }
      if (fwrite(in.data(), sizeof(T), count, out) != count)
      {
        fclose(out);
        unlink(filename.c_str());
        throw std::runtime_error("Error writing '" + filename + "'!");
      }
      fclose(out);
    }

  } // end namespace
} // end namespace

#endif
