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

/*! \file src/Healpix_2.15a/cxxsupport/cxxutils.h
 *  Various convenience functions used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke \author Reinhard Hell
 */

#ifndef PLANCK_CXXUTILS_H
#define PLANCK_CXXUTILS_H

#include <algorithm>
#include <string>
#include <map>
#include <cmath>
#include "src/Healpix_2.15a/message_error.h"
#include "src/Healpix_2.15a/lsconstants.h"

/*! \defgroup mathutilsgroup Mathematical helper functions */
/*! \{ */

//! Returns \e true if | \a a-b | < \a epsilon * | \a b |, else \e false.
template<typename F> inline bool approx (F a, F b, F epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) < (epsilon*abs(b));
  }

//! Returns \e true if | \a a-b | < \a epsilon, else \e false.
template<typename F> inline bool abs_approx (F a, F b, F epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) < epsilon;
  }

//! Returns the largest integer which is smaller than (or equal to) \a arg.
template<typename I, typename F> inline I ifloor (F arg)
  {
  return (arg>=0) ? I(arg) : I(arg)-1;
  }

//! Returns the integer which is nearest to \a arg.
template<typename I, typename F> inline I nearest (F arg)
  {
  arg += 0.5;
  return (arg>=0) ? I(arg) : I(arg)-1;
  }

//! Returns \a v1+v2 if \a v1<0, \a v1-v2 if \a v1>=v2, else \a v1.
/*! \a v1 can be positive or negative; \a v2 must be positive. */
template<typename T> inline T weak_modulo (T v1, T v2)
  { return (v1>=0) ? ((v1<v2) ? v1 : (v1-v2)) : (v1+v2); }

//! Returns the remainder of the division \a v1/v2.
/*! The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
inline double fmodulo (double v1, double v2)
  {
  using namespace std;
  return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
  }

//! Returns the remainder of the division \a v1/v2.
/*! The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
template<typename I> inline I imodulo (I v1, I v2)
  { return (v1>=0) ? ((v1<v2) ? v1 : (v1%v2)) : ((v1%v2)+v2); }

//! Returns -1 if \a signvalue is negative, else +1.
template<typename T> inline T sign (const T& signvalue)
  { return (signvalue>=0) ? 1 : -1; }

//! Returns the integer \a n, which fulfills \a n*n<=arg<(n+1)*(n+1).
template<typename I> inline unsigned int isqrt (I arg)
  {
	  using namespace std;
	  if (sizeof(I)<=4)
		return unsigned (sqrt(arg+0.5));
	  else
	  {
		  long double arg2 = arg;
		  return unsigned (sqrt(arg2+0.5));
	  }
  }

//! Returns the largest integer \a n that fulfills \a 2^n<=arg.
template<typename I> inline unsigned int ilog2 (I arg)
  {
  unsigned int res=0;
  while (arg > 0x0000FFFF) { res+=16; arg>>=16; }
  if (arg > 0x000000FF) { res|=8; arg>>=8; }
  if (arg > 0x0000000F) { res|=4; arg>>=4; }
  if (arg > 0x00000003) { res|=2; arg>>=2; }
  if (arg > 0x00000001) { res|=1; }
  return res;
  }

//! Returns \a atan2(y,x) if \a x!=0 or \a y!=0; else returns 0.
inline double safe_atan2 (double y, double x)
  {
  using namespace std;
  return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x);
  }

//! Returns an index to the left of two interpolation values.
/*! \a begin points to an array containing a sequence of values
    sorted in ascending order. The length of the array is \a len.
    If \a val is lower than the first element, 0 is returned.
    If \a val is higher than the last element, \a len-2
    is returned. Else, the index of the largest element smaller
    than \a val is returned. */
template<typename T> inline int interpol_left
  (const T *begin, int len, const T &val)
  {
  const T *end = begin+len;
  const T *iter = std::lower_bound (begin, end, val);
  if (iter==begin) return 0;
  if (iter==end) return len-2;
  return (iter-begin)-1;
  }

//! Returns an index to the nearest interpolation value.
/*! \a begin points to an array containing a sequence of values
    sorted in ascending order. The length of the array is \a len.
    If \a val is lower than the first element, 0 is returned.
    If \a val is higher than the last element, \a len-1 is returned.
    Else, the index of the nearest element within the sequence of
    values is returned. */
template<typename T> inline int interpol_nearest
  (const T *begin, int len, const T &val)
  {
  int left = interpol_left(begin, len, val);
  T delleft = val-(*(begin+left));
  T delright = (*(begin+left+1))-val;
  if (delright<0) return left+1;
  return (delright<delleft) ? (left+1) : left;
  }

/*! \} */

/*! \defgroup fileutilsgroup File-handling helper functions */
/*! \{ */

//! If the file \a filename is present, return \p true, else \p false.
bool file_present (const std::string &filename);

//! Removes the file \a filename
void remove_file (const std::string &filename);

/*! \} */

/*! \defgroup assertgroup Assertions */
/*! \{ */

//! Throws a Message_error containing \a msg if \a testval is false.
inline void planck_assert (bool testval, const std::string &msg)
  {
  if (testval) return;
  throw Message_error ("Assertion failed: "+msg);
  }
//! Throws a Message_error containing \a msg if \a testval is false.
inline void planck_assert (bool testval, const char *msg)
  {
  if (testval) return;
  throw Message_error ("Assertion failed: "+std::string(msg));
  }

//! Checks the presence of the file \a filename.
/*! If the file is not present, a Message_error is thrown. */
void assert_present (const std::string &filename);

//! Checks the absence of the file \a filename.
/*! If the file is present, a Message_error is thrown. */
void assert_not_present (const std::string &filename);

/*! \} */

/*! \defgroup stringutilsgroup String handling helper functions */
/*! \{ */

//! Returns the string \a orig without leading and trailing whitespace.
std::string trim (const std::string &orig);

//! Returns a string containing the text representation of \a x.
/*! Care is taken that no information is lost in the conversion. */
template<typename T> std::string dataToString(const T &x);
template<> std::string dataToString (const bool &x);
template<> std::string dataToString (const std::string &x);
template<> std::string dataToString (const float &x);
template<> std::string dataToString (const double &x);

/*! Returns a string containing the text representation of \a x, padded
    with leading zeroes to \a width characters. */
std::string intToString(int x, int width);

//! Reads a value of a given datatype from a string
template<typename T> void stringToData (const std::string &x, T &value);
template<> void stringToData (const std::string &x, std::string &value);
template<> void stringToData (const std::string &x, bool &value);

//! Reads a value of a given datatype from a string
template<typename T> inline T stringToData (const std::string &x)
  { T result; stringToData(x,result); return result; }

//! Parses the file \a filename and returns the key/value pairs in \a dict.
void parse_file (const std::string &filename,
  std::map<std::string,std::string> &dict);

//! Case-insensitive string comparison
/*! Returns \a true, if \a a and \a b differ only in capitalisation,
    else \a false. */
bool equal_nocase (const std::string &a, const std::string &b);

//! Returns lowercase version of \a input.
std::string tolower(const std::string &input);

/*! \} */

//! Indicates progress by printing the percentage of \a now/total.
/*! A message is only printed if it has changed since \a now-1/total.
    The output is followed by a carriage return, not a newline. */
void announce_progress (int now, int total);
//! Indicates progress by printing the percentage of \a now/total.
/*! A message is only printed if it has changed since \a last/total.
    The output is followed by a carriage return, not a newline. */
void announce_progress (double now, double last, double total);
/*! This function should be called after a sequence of announce_progress()
    calls has finished. */
void end_announce_progress ();

//! Prints a banner containing \a name. Useful for displaying program names.
void announce (const std::string &name);

/*! Prints a banner containing \a name and checks if \a argc==argc_expected.
    If not, a usage description is given and the program is terminated. */
void module_startup (const std::string &name, int argc, const char **argv,
  int argc_expected, const std::string &argv_expected);

//! Returns an appropriate FITS repetition count for a map with \a npix pixels.
inline unsigned int healpix_repcount (int npix)
  {
  if (npix<1024) return 1;
  else if ((npix%1024)==0) return 1024;
  else return isqrt (npix/12);
  }
#endif
