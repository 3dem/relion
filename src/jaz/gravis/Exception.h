#ifndef __LIBGRAVIS_EXCEPTION_H__
#define __LIBGRAVIS_EXCEPTION_H__
/******************************************************************************
**	  Title: Exception.h
**  Description: Base class for exceptions in libgravis.
**
**	 Author: Jean-Sebastien Pierrard, 2005
**		 Computer Science Department, University Basel (CH)
**
******************************************************************************/

#define GRAVIS_CHECK(condition, r) if (!(condition)) GRAVIS_THROW3(gravis::Exception, "Assertion failed", #condition)
#define GRAVIS_THROW(e)	       throw e(std::string(__FILE__),__LINE__)
#define GRAVIS_THROW2(e,r)     throw e(std::string(__FILE__),__LINE__,(r))
#define GRAVIS_THROW3(e,r,arg) throw e(std::string(__FILE__),__LINE__,(r),(arg))

#include <string>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>

namespace gravis
{

  class Exception : public std::runtime_error
  {

    public:

      Exception (const std::string& src, const int line,
                 const std::string& dtl = "", const std::string& arg = "")
        : std::runtime_error( std::string("gravis exception: ") + src + ", " + dtl + " (" + arg + ")" ),
          _source(src), _detail(dtl), _arg(arg), _line(line) {}

      Exception(const Exception& e) : std::runtime_error( e.what() ), _source(e._source), _detail(e._detail), _arg(e._arg), _line(e._line) { }

      virtual ~Exception() throw() {}

      virtual const char* getClassName () const
      {
        return "Exception";
      }

      const char* detail() const
      {
        return _detail.c_str();
      }
      const char* argument() const
      {
        return _arg.c_str();
      }
      const char* source() const
      {
        return _source.c_str();
      }
      const int& line() const
      {
        return _line;
      }
      bool hasDetail() const
      {
        return _detail.length() > 0;
      }
      bool hasArgument() const
      {
        return _arg.length() > 0;
      }

      /* do not need that since ctor of runtime_error took the message already
      virtual const char* what() const throw()
      {
        std::string strError( _source + " " + _detail + " " + _arg );

        char* pLostMem = new char[ strError.size() + 1 ];
        for( size_t i = 0; i < strError.size(); i++ )
          pLostMem[i] = strError[i];

        pLostMem[ strError.size() ] = '\0';
        return pLostMem;
      }*/

    protected:

      std::string _source;
      std::string _detail;
      std::string _arg;
      int _line;
  };

  inline
  std::ostream& operator << (std::ostream& os, const Exception& ex)
  {
    os << ex.getClassName() << " in  " << ex.source() << ", line " << ex.line();
    if (ex.hasDetail()) os << ": " << ex.detail();
    if (ex.hasArgument()) os << " (" << ex.argument() << ")";
    return os;
  }

} /* Close Namespace "gravis" */

#endif
