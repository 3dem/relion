/******************************************************************************
**        Title: gravis/io/array.h
**  Description: Implements reader/writer for different array file formats.
**
**       Author: Sandro Schoenborn
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef __GRAVIS_IO_ARRAY__
#define __GRAVIS_IO_ARRAY__

#include "array/a.h"
#include "array/ArrayStreamIO.h"
#include "array/le.h"

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace gravis
{
  namespace io
  {

    /**
     * The main functionality of this class is to forward the commands to
     * ArrayA, ArrayLittleEndian, ArrayStreamIO based on the file contents.
     *
     **/
    class Array
    {

      private:
        static
        inline
        bool has_ending(const std::string& filename, const std::string& ending)
        {
          return boost::algorithm::ends_with( filename, ending );
        }

        static
        inline
        bool is_readable_format( std::string const& filename )
        {
          if ( has_ending( filename, ".gz" ) )
          { // compressed format -- so far only readable in this format
            return true;
          }
          else
          { // uncompressed format
            std::ifstream ifIn( filename.c_str() );
            std::string strLine;
            std::getline( ifIn, strLine );

            if ( strLine.empty() )
              return false;

            boost::algorithm::erase_all( strLine, " " );
            size_t nP1 = strLine.find("[0]"); // first index is alway 0
            if ( nP1 == 0 )
              return true;
            else
              return false;
          }
        }

        static
        inline
        bool is_header_format( std::string const& filename )
        {
          std::ifstream ifIn( filename.c_str(), std::ios_base::binary );

          std::string strMagic = ArrayA::magic();
          std::vector<char> vCheck( strMagic.size() );

          if (  ifIn.read( &vCheck[0], vCheck.size() ) )
          {
            std::string strCheck( vCheck.begin(), vCheck.end() );

            if ( strCheck == strMagic )
              return true;
            else
              return false;
          }
          else
            return false;
        }

      public:
         /**
         * Load vector from a file. The filetype is determined automatically
         **/
        template<typename T>
        static
        inline
        std::vector<T> load( const std::string& filename )
        {
            std::vector<T> vec;
            load( vec, filename );
            return std::move(vec);
        }

        /**
         * Load vector from a file. The filetype is determined automatically
         **/
        template<typename T>
        static
        inline
        void load( std::vector<T>& out, const std::string& filename, size_t count = 0 )
        {
          if ( is_readable_format( filename ) )
            ArrayStreamIO::load( out, filename );
          else if ( is_header_format( filename ) )
            ArrayA::load( out, filename );
          else
          {
            tArray<T> tout;
            ArrayLittleEndian::load( tout, filename, count );
            out = std::vector<T>( tout );
          }
        }

        /**
         * save array to a file, standard format is ArrayStreamIO (readable format)
         **/
        template<typename T>
        static
        inline
        void save(const std::string& filename, const std::vector<T>& array )
        {
          ArrayStreamIO::save( filename, array );
        }

        /**
         * Load vector from a file. The filetype is determined automatically
         **/
        template<typename T>
        static
        inline
        void load( gravis::tArray<T>& out, const std::string& filename, size_t count = 0 )
        {
          if ( is_readable_format( filename ) )
            ArrayStreamIO::load( out, filename );
          else if ( is_header_format( filename ) )
            ArrayA::load( out, filename );
          else
            ArrayLittleEndian::load( out, filename, count );
        }

        /**
         * save array to a file, standard format is ArrayStreamIO (readable format)

         **/
        template<typename T>
        static
        inline
        void save(const std::string& filename, const gravis::tArray<T>& array )
        {
          ArrayStreamIO::save( filename, array );
        }
    };

  }
}

#endif
