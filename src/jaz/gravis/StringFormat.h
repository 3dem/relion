#ifndef __LIBGRAVIS_STRING_FORMAT_H__
#define __LIBGRAVIS_STRING_FORMAT_H__

#include <sstream>
#include <iomanip>
#include <stdint.h>

namespace gravis
{
  /**
   * Usage
   *  std::string name = StringFormat("Dies sind ")(12)(" Zahlen.");
   *  std::string name = StringFormat("Dies sind ")(12, 4, '0')(" Zahlen.");
   **/
  class StringFormat
  {

    private:
      std::stringstream s;

    public:
      StringFormat(const StringFormat& start) : s()
      {
        s << start.string();
      };

// Takanori removed this because this is broken and not used.
/*    const char*   c_str()  const
      {
        return s.str().c_str();
      }
*/
      std::string   string() const
      {
        return s.str();
      }
      operator std::string() const
      {
        return s.str();
      }

      bool operator==(const StringFormat& o) const
      {
        return o.string()==string();
      }
      bool operator!=(const StringFormat& o) const
      {
        return o.string()!=string();
      }
      bool operator==(const std::string& o) const
      {
        return o==string();
      }
      bool operator!=(const std::string& o) const
      {
        return o!=string();
      }

      StringFormat() : s() {  }

      template <class T>
      explicit StringFormat(const T&                      e) : s()
      {
        s << e;
      }

      template <class T>
      StringFormat(const T& e, std::streamsize w) : s()
      {
        s << std::setw(w) << e;
      }

      template <class T>
      StringFormat(const T&   e, std::streamsize w, char fill) : s()
      {
        s << std::setw(w) << std::setfill(fill) << e;
      }

      template <class T>
      inline StringFormat& operator()(const T& e)
      {
        s << e;
        return *this;
      }

      template <class T>
      inline StringFormat& operator()(const T& e, int w)
      {
        s << std::setw(w) << e;
        return *this;
      }

      template <class T>
      inline StringFormat& operator()(const T& e, int w, char fill)
      {
        s << std::setw(w) << std::setfill(fill) << e;
        return *this;
      }
  };

}



inline
std::ostream& operator<< (std::ostream& os, const gravis::StringFormat& arg)
{
  os << arg.string();
  return os;
}

#endif
