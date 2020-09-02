/******************************************************************************
 **        Title: io/array/a.h
 **  Description: Input/Output for tArrays.
 **
 **       Author: Brian Amberg, 2006
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#ifndef __GRAVIS_IO_ARRAY_A__
#define __GRAVIS_IO_ARRAY_A__

#include "../../tArray.h"
#include "../../Exception.h"
#include <fstream>
#include <string>
#include <cstring>
#include <stdint.h>
#include "../../Tuple.h"
#include "../../t2Vector.h"
#include "../../t3Vector.h"
#include "../../t4Vector.h"
#include "../../StringFormat.h"

namespace gravis
{

  namespace io
  {

    /** This should be in ArrayA in the private section, but that is not allowed **/
    namespace ArrayA_private
    {

      template <class T> struct TName
      {
        static std::string name()
        {
          return "UNKNOWN";
        }
      };
      template <>        struct TName<float>
      {
        static std::string name()
        {
          return "float";
        }
      };
      template <>        struct TName<double>
      {
        static std::string name()
        {
          return "double";
        }
      };
      template <>        struct TName<char>
      {
        static std::string name()
        {
          return "char";
        }
      };
      template <>        struct TName<int8_t>
      {
        static std::string name()
        {
          return "int8";
        }
      };
      template <>        struct TName<int16_t>
      {
        static std::string name()
        {
          return "int16";
        }
      };
      template <>        struct TName<int32_t>
      {
        static std::string name()
        {
          return "int32";
        }
      };
      template <>        struct TName<int64_t>
      {
        static std::string name()
        {
          return "int64";
        }
      };
      template <>        struct TName<uint8_t>
      {
        static std::string name()
        {
          return "uint8";
        }
      };
      template <>        struct TName<uint16_t>
      {
        static std::string name()
        {
          return "uint16";
        }
      };
      template <>        struct TName<uint32_t>
      {
        static std::string name()
        {
          return "uint32";
        }
      };
      template <>        struct TName<uint64_t>
      {
        static std::string name()
        {
          return "uint64";
        }
      };
      template <>        struct TName<Tuple2>
      {
        static std::string name()
        {
          return "uint32";
        }
      };
      template <>        struct TName<Tuple3>
      {
        static std::string name()
        {
          return "uint64";
        }
      };
      template <class T> struct TName<t2Vector<T> >
      {
        static std::string name()
        {
          return std::string("vector2:") + TName<T>::name();
        }
      };
      template <class T> struct TName<t3Vector<T> >
      {
        static std::string name()
        {
          return std::string("vector3:") + TName<T>::name();
        }
      };
      template <class T> struct TName<t4Vector<T> >
      {
        static std::string name()
        {
          return std::string("vector4:") + TName<T>::name();
        }
      };


    }

    /**
     * Simple Array format with a header describing the array type
     **/
    class ArrayA
    {

      private:
        static const uint8_t version = 1;

      public:
        static const std::string magic()
        {
          return "GVS:ARR";
        };

      public:

        struct ArrayHeader
        {
          char magic[8];                  // 0
          char type[16];                  // 8
          uint16_t type_size;             // 24
          uint32_t length;                // 26
          uint8_t version;                // 30
          char reserved1[1];               // 31
          char reserved2[32];              // 32
        };

        /**
         * Load an array from an array file. Type checking is done for a subset
         * of supported types. More types can be added by changing the private
         * part of this class. Unknown types are saved and loaded, but not type checked.
         **/
        template<class T>
        static
        bool is_a(const std::string& filename)
        {
          std::ifstream is(filename.c_str(),  std::ios_base::binary);
          ArrayHeader h;
          is.read((char*)(&h), sizeof(h));
          if (!is.good())
            GRAVIS_THROW3(::gravis::Exception, "Could not read file", filename);
          if (std::string(h.magic) != magic())
            GRAVIS_THROW3(::gravis::Exception, "Not a gravis array file", filename);
          if (h.version > version)
            GRAVIS_THROW3(::gravis::Exception, "Can't read this gravis array version", filename);
          if (sizeof(T) != h.type_size)
            return false;
          if (ArrayA_private::TName<T>::name() != h.type)
            return false;
          return true;
        }

        /**
         * Load an array from an array file. Type checking is done for a subset
         * of supported types. More types can be added by changing the private
         * part of this class. Unknown types are saved and loaded, but not type checked.
         **/
        template<class T>
        static
        void load(std::vector<T> &out, const std::string& filename)
        {
          std::ifstream is(filename.c_str(),  std::ios_base::binary);
          ArrayHeader h;
          is.read((char*)(&h), sizeof(h));
          if (!is.good())
            GRAVIS_THROW3(::gravis::Exception, "Could not read file", filename);
          if (std::string(h.magic) != magic())
            GRAVIS_THROW3(::gravis::Exception, "Not a gravis array file", filename);
          if (h.version > version)
            GRAVIS_THROW3(::gravis::Exception, "Can't read this gravis array version", filename);
          if ((ArrayA_private::TName<T>::name() != h.type) || (sizeof(T) != h.type_size))
            GRAVIS_THROW3(::gravis::Exception, ::gravis::StringFormat("Wrong type in array. Expected ")(ArrayA_private::TName<T>::name())(" of size ")(sizeof(T))(" but got ")(h.type)(" of size ")(h.type_size), filename);
          out.resize(h.length);
          is.read((char*)(&out[0]), sizeof(out[0]) * h.length);
        }

        /**
         * Load an array from an array file. Type checking is done for a subset
         * of supported types. More types can be added by changing the private
         * part of this class. Unknown types are saved and loaded, but not type checked.
         **/
        template<class T>
        static
        void load(tArray<T> &out, const std::string& filename)
        {
          std::ifstream is(filename.c_str(),  std::ios_base::binary);
          ArrayHeader h;
          is.read((char*)(&h), sizeof(h));
          if (!is.good())
            GRAVIS_THROW3(::gravis::Exception, "Could not read file", filename);
          if (std::string(h.magic) != magic())
            GRAVIS_THROW3(::gravis::Exception, "Not a gravis array file", filename);
          if (h.version > version)
            GRAVIS_THROW3(::gravis::Exception, "Can't read this gravis array version", filename);
          if ((ArrayA_private::TName<T>::name() != h.type) || (sizeof(T) != h.type_size))
            GRAVIS_THROW3(::gravis::Exception, ::gravis::StringFormat("Wrong type in array. Expected ")(ArrayA_private::TName<T>::name())(" of size ")(sizeof(T))(" but got ")(h.type)(" of size ")(h.type_size), filename);
          out.resize(h.length);
          is.read((char*)(&out[0]), sizeof(out[0]) * h.length);
        }

        /**
         * Save an array to an array file. Type checking is done for a subset
         * of supported types. More types can be added by changing the private
         * part of this class
         **/
        template<class T>
        static
        void save(const std::string& filename, const tArray<T> &in)
        {
          std::ofstream os(filename.c_str(),  std::ios_base::binary);
          ArrayHeader h;
          std::memset((char*)(&h), 0x0, sizeof(h));
          std::memcpy(h.magic, magic().data(), magic().size());
          std::memcpy(h.type, ArrayA_private::TName<T>::name().data(), ArrayA_private::TName<T>::name().size());
          h.type_size = sizeof(T);
          h.length = in.size();
          h.version = version;
          os.write((char*)(&h), sizeof(h));
          os.write((char*)(&in[0]), sizeof(T) * in.size());
        }


        /**
         * Save an array to an array file. Type checking is done for a subset
         * of supported types. More types can be added by changing the private
         * part of this class
         **/
        template<class T>
        static
        void save(const std::string& filename, const std::vector<T> &in)
        {
          std::ofstream os(filename.c_str(),  std::ios_base::binary);
          ArrayHeader h;
          std::memset((char*)(&h), 0x0, sizeof(h));
          std::memcpy(h.magic, magic().data(), magic().size());
          std::memcpy(h.type, ArrayA_private::TName<T>::name().data(), ArrayA_private::TName<T>::name().size());
          h.type_size = sizeof(T);
          h.length = in.size();
          h.version = version;
          os.write((char*)(&h), sizeof(h));
          os.write((char*)(&in[0]), sizeof(T) * in.size());
        }

    };

  }
}

#endif
