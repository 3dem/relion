/*============================================================================*/
/**
 *         @file ClearTextIO.h
 *
 *        @brief Gravis Arrays with >> and << reader and writer
 *
 *         @date 06/30/2011 02:18:07 PM
 *      @authors Sandro Schoenborn (ses)\n
 *               sandro.schoenborn@unibas.ch\n
 *               University of Basel, Switzerland
 */
/*============================================================================*/


#ifndef  ARRAYSTREAMIO_INC
#define  ARRAYSTREAMIO_INC

#include "../../tArray.h"

#include <stdexcept>

#include <iostream>
#include <sstream>
#include <fstream>
#include <deque>
#include <algorithm>
#include <limits>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace gravis
{
  /** Stream output operator for vectors, each index is explicitely written */
  template<class T>
  std::ostream& operator<< ( std::ostream& ostr, std::vector<T> const& in )
  {
    //ostr << std::scientific;
    ostr.precision(17); // 16 might be just enough

    for ( size_t i = 0; i < in.size(); i++ )
      ostr << "[" << i << "]=" << in[i] << "\n";

    return ostr;
  }

  /** Stream input operator for vectors, each index is explicitely expected, format: "[i]=value\n",
    * if \param in already has a size, it is only filled, not enlarged!
    */
  template<class T>
  std::istream& operator>> ( std::istream& istr, std::vector<T>& in )
  {
    std::string strLine;
    size_t nArrPos = 0;
    size_t nMaxSize = std::numeric_limits<size_t>::max();

    if ( !in.empty() )
    {
      nMaxSize = in.size();
      nArrPos = in.size();
    }

    while ( std::getline( istr, strLine ) && nArrPos <= nMaxSize )
    {
      size_t nPos = 0;
      T tCont;

      // Do not parse empty lines
      if ( strLine.empty() )
        continue; 

      // expect format: [index]=value
      size_t nEq = strLine.find('=');
      size_t nB1 = strLine.find('[');
      size_t nB2 = strLine.find(']');

      if ( nEq == std::string::npos || nB1 == std::string::npos || nB2 == std::string::npos
           ||  nEq < nB1 || nEq < nB2 || nB2 < nB1 )
        throw std::runtime_error("ArrayStreamIO: format is invalid! expect: \"[index]=value\" per line");

      std::string strIndex = strLine.substr( nB1+1, nB2-nB1-1 );
      boost::algorithm::trim( strIndex );

      std::string strContent = strLine.substr( nEq+1 );

      // Index number
      std::stringstream( strIndex ) >> nPos;
      // Content
      std::stringstream( strContent ) >> tCont;

      if ( nPos == nArrPos )
      {
        in.push_back( tCont ); // deque push_back should not be too costly
      }
      else if ( nPos < nArrPos )
      {
        in[nPos] = tCont;
      }
      else
      {
        //std::cout << "Larger index, resizing to " << nPos << std::endl;
        in.resize( nPos+1 );
        in[nPos] = tCont;
      }
      nArrPos = in.size();
    }

    if ( in.size() > nMaxSize && nMaxSize < std::numeric_limits<size_t>::max() )
      in.resize( nMaxSize );

    return istr;
  }

  namespace io
  {
    /** Provides vector (and tArray) IO via streams, allows for gzip compression */
    class ArrayStreamIO
    {
      public:
        /** Load an array (vector) from a file, using stream operator>>, able to decompress gzip files (filename ends with gz)*/
        template<class T>
        static void load( std::vector<T>& out, std::string const& filename);

        /** Save an array (vector) to a file, using stream operator<<, able to compress gzip files (filename ends with gz)*/
        template<class T>
        static void save( std::string const& filename, std::vector<T> const& in, bool compression = false );

        /** Load an array (tArray) from a file, using stream operator>>, able to decompress gzip files (filename ends with gz)*/
        template<class T>
        static void load( tArray<T>& out, std::string const& filename);

        /** Save an array (tArray) to a file, using stream operator<<, able to compress gzip files (filename ends with gz)*/
        template<class T>
        static void save( std::string const& filename, tArray<T> const& in, bool compression = false );
    };

    // ------ Implementation ------
    template<class T>
    void ArrayStreamIO::load(std::vector<T>& out, std::string const& filename )
    {
      if ( boost::algorithm::ends_with( filename, "gz" ) )
      {
        // compressed file read with gzip filter
        std::ifstream ifInput( filename.c_str(), std::ios_base::in | std::ios_base::binary );
        if (!ifInput)
          throw std::runtime_error("Cannot open '" + filename + "'!");
        boost::iostreams::filtering_streambuf<boost::iostreams::input> bufIn;
        bufIn.push( boost::iostreams::gzip_decompressor() );
        bufIn.push( ifInput );
        std::stringstream ssData;
        boost::iostreams::copy( bufIn, ssData );
        ssData >> out;
      }
      else
      { // normal readable file
        std::ifstream ifInput( filename.c_str(), std::ios_base::in );
        if (!ifInput)
          throw std::runtime_error("Cannot open '" + filename + "'!");
        ifInput >> out;
      }
    }

    template<class T>
    void ArrayStreamIO::save(std::string const& filename, std::vector<T> const& in, bool compression )
    {
      //      if (!count) count = in.size();
      if ( compression || boost::algorithm::ends_with( filename, "gz" ) )
      {
        std::string strFile( filename );
        if ( !boost::algorithm::ends_with( strFile, "gz" ) )
          strFile += ".gz";

        std::ofstream ofOut( strFile.c_str(), std::ios_base::out | std::ios_base::binary );
        if ( !ofOut )
          throw std::runtime_error( "Could not open file for writing: " + strFile );

        boost::iostreams::filtering_streambuf<boost::iostreams::input> bufOut;
        bufOut.push( boost::iostreams::gzip_compressor() ); // write a gzip compressed file

        std::stringstream ssData;
        ssData << in; // via stream operator
        bufOut.push( ssData );

        boost::iostreams::copy( bufOut, ofOut );
      }
      else
      {
        std::ofstream ofOut( filename.c_str(), std::ios_base::out );
        if ( !ofOut )
          throw std::runtime_error( "Could not open file for writing: " + filename );

        ofOut << in; // via stream operator
      }
    }

    // tArray compatibility stuff
    template<class T>
    void ArrayStreamIO::load(tArray<T>& out, std::string const& filename )
    {
      std::ifstream ifInput( filename.c_str(), std::ios_base::in );

      if (!ifInput)
      {
        throw std::runtime_error("Cannot open '" + filename + "'!");
      }

      std::vector<T> vOut( out.size() );
      load( vOut, filename );
      out = tArray<T>( vOut );
    }

    template<class T>
    void ArrayStreamIO::save(std::string const& filename, tArray<T> const& in, bool compression )
    {
      std::vector<T> vArr( in );
      save( filename, vArr, compression );
    }
  }
}



#endif   // ----- #ifndef ARRAYSTREAMIO_INC  -----
