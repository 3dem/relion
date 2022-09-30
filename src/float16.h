/***************************************************************************
 *
 * Author: "Takanori Nakane"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef FLOAT16_H
#define FLOAT16_H

/*
   BIT LAYOUT FOR IEEE 754 FLOATING POINT NUMBERS

   Reference: https://en.wikipedia.org/wiki/Single-precision_floating-point_format
              https://en.wikipedia.org/wiki/Half-precision_floating-point_format

   - float32:
      1 + 8 + 23 bits, exp offset = 127 (0..254 correspond to -126 to 127; 255 is INF/NAN)
      seee eeee | efff ffff | ffff ffff | ffff ffff
      sign mask = 0x80000000, exp mask = 0x7f800000, frac mask = 0x007fffff
   - float16:
      1 + 5 + 10 bits, exp offset = 15 (0..30 correspond to -14 to 15; 31 is INF/NAN)
      seee eeff | ffff ffff
      sign mask = 0x8000, exp mask = 0x7c00, frac mask = 0x03ff

   SPECIAL NUMBERS

   - (Signed) zero: sign = any, exp = 0, frac = 0
   - Subnormal number: sign = any, exp = 0, frac != 0
   - Infinity:  sign = any, exp = 255 (or 31), frac = 0
   - NaN:  sign = any, exp = 255 (or 31), frac != 0
*/

// We cannot simply cast pointers; it will break the strict aliasing rule.
//
// Strictly speaking, type-punning via union is allowed in C11 but
// appears to be an undefined behavior in C++11.
// https://stackoverflow.com/questions/11373203/accessing-inactive-union-member-and-undefined-behavior
// However, at least GCC and ICC support this. (Am I correct??)
// http://gcc.gnu.org/onlinedocs/gcc-4.7.1/gcc/Structures-unions-enumerations-and-bit_002dfields-implementation.html#Structures-unions-enumerations-and-bit_002dfields-implementation
// https://www.intel.com/content/www/us/en/developer/articles/technical/pointer-aliasing-and-vectorization.html
typedef union
{
	float f;
	unsigned int i;
} float32;

typedef unsigned short float16;

/*
  Conversion to and from float32.
 
  I learned the algorithm from NumPy's implementation (BSD license)
  https://github.com/numpy/numpy/blob/master/numpy/core/src/npymath/halffloat.c
  but wrote this myself from scratch.

  My code below is more naive (but less cryptic!) and hopefully compilers can optimize it.
  An important difference is that we don't want to generate INFs. Overflows and underflows
  are truncated to the maximum and minimum normal numbers float16 can represent.

  At the moment subnormal numbers are not implemented; they are truncated to signed zeros.
*/

inline float16 float2half(float f32)
{
	float32 f;
	f.f = f32;

	unsigned int sign = (f.i & 0x80000000u) >> 31; // 1 bit
	unsigned int exponent = (f.i & 0x7f800000u) >> 23; // 10 bits
	unsigned int fractional = (f.i & 0x007fffffu); // 23 bits

	unsigned short ret = sign << 15; // first make a signed zero

	if (exponent == 0)
	{
		// Do nothing. Subnormal numbers will be signed zero.
		return ret;
	}
	else if (exponent == 255) // Inf, -Inf, NaN 
	{
		ret |= 31 << 10; // exponent is 31
		ret |= fractional >> 13; // fractional is truncated from 23 bits to 10 bits
		return ret;
	}

	fractional += 1 << 12; // add 1 to 13th bit to round.
	if (fractional & (1 << 23)) // carry up
		exponent++;

	if (exponent > 127 + 15) // Overflow: don't create INF but truncate to MAX.
	{
		ret |= 30 << 10; // maximum exponent 30 (= +15)
		ret |= 0x03ffu; // 10 bits of 1s
		return ret;
	}
	else if (exponent < 127 - 14) // Underflow
	{
		return ret; // TODO: generate subnormali numbers instead of returning a signed zero
	}
	else
	{
		ret |= ((exponent + 15 - 127) & 0x1f) << 10;
		ret |= fractional >> 13; // fractional is truncated from 23 bits to 10 bits
		return ret;
	}
}

inline float half2float(float16 f16)
{
	unsigned int sign = (f16 & 0x8000u) >> 15; // 1 bit
	unsigned int exponent = (f16 & 0x7c00u) >> 10; // 5 bits
	unsigned int fractional = (f16 & 0x03ffu); // 10 bits

	unsigned int ret = sign << 31; // first make a signed zero

	if (exponent == 0)
	{
		if (fractional != 0) {
			// TODO: convert float16 subnormal numbers to float32 normal numbers
		}
		// else signed zero: nothing to do
	}
	else if (exponent == 31) // Inf, -Inf, NaN
	{
		ret |= 255 << 23; // exponent is 255
		ret |= fractional << 13; // keep fractional by expanding from 10 bits to 23 bits
	}
	else // normal numbers
	{
		ret |= (exponent + 127 - 15) << 23; // shift the offset
		ret |= fractional << 13; // keep fractional by expanding from 10 bits to 23 bits
	}

	float32 f;
	f.i = ret;

	return f.f;
}

#endif
