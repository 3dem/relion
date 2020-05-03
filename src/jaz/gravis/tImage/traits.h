/******************************************************************************
**        Title: tImage/traits.h
**  Description: tImage traits neccessary for nice implementation of stuff
**               like convolution
**
**       Author: Brian Amberg
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef _TIMAGE_TRAITS_H_
#define _TIMAGE_TRAITS_H_

// TODO: This centralizes traits definition, but includes everything, which is
// bad. It would be better to use predeclarations of these classes here, but
// that does not work with templates. Another possibility would be to move the
// trait definitions into the include files of the respective datatypes.

#include "../t2Vector.h"
#include "../t3Vector.h"
#include "../t4Vector.h"
#include "../t2Matrix.h"
#include "../t3Matrix.h"
#include "../t4Matrix.h"
#include "../tMatrix.h"
#include "../tRGB.h"
#include "../tBGR.h"
#include "../tYCbCr.h"
#include "../tRGBA.h"
#include "../tRGB_A.h"
#include "../tGray_A.h"
#include "../tLab.h"

namespace gravis
{

  template <class T>
  struct tImageTraits
  {
      typedef float Scalar_t;
      typedef float Float_t;
      typedef float Pixel_t;
  };

#define DEFINE_TRAIT( aPixel_t, aScalar_t, aFloat_t) \
  template <>                                        \
  struct tImageTraits< aPixel_t > {                  \
    typedef aScalar_t    Scalar_t;                   \
    typedef aFloat_t     Float_t;                    \
    typedef aPixel_t     Pixel_t;                    \
  	static unsigned int components(){          \
			return sizeof(aPixel_t)/sizeof(Scalar_t);      \
		}                                                \
  }

#define DEFINE_ALL_COMPOUND_TRAITS( aScalar_t, aFloat_t) \
  DEFINE_TRAIT( aScalar_t,             aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tRGB< aScalar_t >,     aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tBGR< aScalar_t >,     aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tRGBA< aScalar_t >,    aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tRGB_A< aScalar_t >,   aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tGray_A< aScalar_t >,  aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tYCbCr< aScalar_t >,   aScalar_t, aFloat_t); \
  DEFINE_TRAIT( tLab< aScalar_t >,     aScalar_t, aFloat_t); \
  DEFINE_TRAIT( t2Vector< aScalar_t >, aScalar_t, aFloat_t); \
  DEFINE_TRAIT( t3Vector< aScalar_t >, aScalar_t, aFloat_t); \
  DEFINE_TRAIT( t4Vector< aScalar_t >, aScalar_t, aFloat_t); \
  DEFINE_TRAIT( t2Matrix< aScalar_t >, aScalar_t, aFloat_t); \
  DEFINE_TRAIT( t3Matrix< aScalar_t >, aScalar_t, aFloat_t); \
  DEFINE_TRAIT( t4Matrix< aScalar_t >, aScalar_t, aFloat_t)


  DEFINE_ALL_COMPOUND_TRAITS( char,                 float );
  DEFINE_ALL_COMPOUND_TRAITS( unsigned char,        float );
  DEFINE_ALL_COMPOUND_TRAITS( signed char,          float );
  DEFINE_ALL_COMPOUND_TRAITS( unsigned int,         double );
  DEFINE_ALL_COMPOUND_TRAITS( signed int,           double );
  DEFINE_ALL_COMPOUND_TRAITS( unsigned short int,   double );
  DEFINE_ALL_COMPOUND_TRAITS( signed short int,     double );
  DEFINE_ALL_COMPOUND_TRAITS( signed long int,      double );
  DEFINE_ALL_COMPOUND_TRAITS( unsigned long int,    double );
  DEFINE_ALL_COMPOUND_TRAITS( float,                float );
  DEFINE_ALL_COMPOUND_TRAITS( double,               double );
  DEFINE_ALL_COMPOUND_TRAITS( long double,          long double );
  DEFINE_ALL_COMPOUND_TRAITS( bool,                 double );

}

#endif
