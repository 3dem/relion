/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Class for error reporting
 *
 *  Copyright (C) 2003, 2004 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_MESSAGE_ERROR_H
#define PLANCK_MESSAGE_ERROR_H

#include <exception>
#include <iostream>
#include <string>

#if defined (PLANCK_STACKTRACE)
#include <execinfo.h>
#endif

inline void show_stackframe()
  {
#if defined (PLANCK_STACKTRACE)
  void *trace[16];
  int trace_size = backtrace(trace, 16);
  char **messages = backtrace_symbols(trace, trace_size);
  std::cerr << "[bt] Execution path:" << std::endl;
  for (int i=0; i<trace_size; ++i)
    std::cerr << "[bt] " << messages[i] << std::endl;
#endif
  }


class Message_error
  {
  private:
    std::string msg;

  public:
    Message_error()
      : msg (std::string("Unspecified error"))
      { std::cerr<<msg<<std::endl; show_stackframe(); }

    explicit Message_error(const std::string &message)
      : msg (message) { std::cerr<<msg<<std::endl; show_stackframe(); }

    virtual const char* what() const
      { return msg.c_str(); }

    virtual ~Message_error() {}
  };

#if defined (PLANCK_CHECKS)

#define PLANCK_DIAGNOSIS_BEGIN try {
#define PLANCK_DIAGNOSIS_END \
} \
catch (Message_error &e) \
  { std::cerr << "Planck exception: " << e.what() << std::endl; throw; } \
catch (std::exception &e) \
  { std::cerr << "std::exception: " << e.what() << std::endl; throw; } \
catch (...) \
  { std::cerr << "Unknown exception" << std::endl; throw; }

#else

#define PLANCK_DIAGNOSIS_BEGIN
#define PLANCK_DIAGNOSIS_END

#endif

#endif
