#ifndef __LIBGRAVIS_PROGRAM_OPTIONS_H__
#define __LIBGRAVIS_PROGRAM_OPTIONS_H__
/**
 * \file
 * This headers includes the boost program options header, disabling the warnings such that we can still have our programs compile with -Werror
 * Also it defines the MACROS PO_SWITCH and PO_VALUE, which make option definitions much more readable.
 **/

#if defined __GNUC__
#pragma GCC system_header
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif

#include <boost/program_options.hpp>
#include <vector>
#include <string>
#include <sstream>

namespace boost
{
  namespace program_options
  {
    template <class T>
    std::string vec2string(const std::vector<T> &v)
    {
      std::stringstream s;
      for (size_t i=0; i<v.size(); ++i) s << v[i];
      return s.str();
    }
  }
}

/**
 * Shortcut to define a boost program options switch.
 * A switch is a boolean value. Use it like this:
 *
 *         bool show_help = false;
 *
 *         po::options_description desc("Video Clicker Options");
 *         desc.add_options()
 *           PO_SWITCH("help,?", show_help, "produce help message")
 *           ;
 **/
#define PO_SWITCH(key, opt, desc) (key, boost::program_options::bool_switch(&opt)->default_value(opt), desc)

/**
 * Shortcut to define a boost program option with a value (e.g. a string or an integer)
 *
 *         int number = false;
 *
 *         po::options_description desc("Example Options");
 *         desc.add_options()
 *           PO_VALUE("number,n", number, "Set the number of foos to use")
 *           ;
 **/
#define PO_VALUE( key, opt, desc) (key, boost::program_options::value(&opt)->default_value(opt), desc)

/**
 * Shortcut to define a boost program option with a vector value (e.g. a vector of strings)
 *
 *         std::vector<int> vec; vec.push_back(1); vec.push_back(2);
 *
 *         po::options_description desc("Example Options");
 *         desc.add_options()
 *           PO_VECTOR("vec,v",  vec,  "List of foos")
 *           ;
 **/
#define PO_VECTOR(key, opt, desc) (key, boost::program_options::value(&opt)->default_value(opt, boost::program_options::vec2string(opt)), desc)

#if defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

#endif
