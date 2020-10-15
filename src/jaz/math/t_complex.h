/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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
#ifndef T_COMPLEX_H
#define T_COMPLEX_H

#include <iostream>
#include <cmath>
//#include "macros.h"

template<class T>
class tComplex
{
	public:

        tComplex()
        {}

        tComplex(T real, T imag = 0)
        :   real(real),
            imag(imag)
        {}


            T real, imag;


        inline tComplex& operator += (const tComplex& arg)
        {
            real += arg.real;
            imag += arg.imag;

            return *this;
        }
		
		inline tComplex& operator -= (const tComplex& arg)
        {
            real -= arg.real;
            imag -= arg.imag;

            return *this;
        }
		
		inline tComplex operator - () const
        {
            return tComplex<T>(-real, -imag);
        }

        inline tComplex& operator *= (const tComplex& arg)
        {
            T re = real*arg.real - imag*arg.imag;
            T im = real*arg.imag + imag*arg.real;

            real = re;
            imag = im;

            return *this;
        }

        inline tComplex& operator /= (const tComplex& arg)
        {
            T cd = arg.real*arg.real + arg.imag*arg.imag;

            T re = real*arg.real + imag*arg.imag;
            T im = imag*arg.real - real*arg.imag;

            real = re/cd;
            imag = im/cd;

            return *this;
        }
		
		inline bool operator == (const tComplex& arg) const
		{
			return (real == arg.real && imag == arg.imag);
		}
  
		inline bool operator != (const tComplex& arg) const
		{
			return !(*this == arg);
		}
		
		inline operator T() const
		{
			return real;
		}

        inline tComplex conj() const
        {
            return tComplex(real, -imag);
        }

        inline T abs() const
        {
            return sqrt(real*real + imag*imag);
        }

        inline T norm() const
        {
            return real*real + imag*imag;
        }

        inline T arg() const
        {
            return atan2(imag,real);
        }
};

template <class T> inline
tComplex<T> conj(const tComplex<T>& op)
{
    return op.conj();
}

template <class T> inline
T abs(const tComplex<T>& op)
{
    return op.abs();
}

template <class T> inline
T norm(const tComplex<T>& op)
{
    return op.norm();
}

template <class T> inline
T arg(const tComplex<T>& op)
{
    return op.arg();
}

template <class T1, class T2> inline
tComplex<T1> operator + (const tComplex<T1>& z, const tComplex<T2>& w)
{
    return tComplex<T1>(z.real + w.real, z.imag + w.imag);
}

template <class T1, class T2> inline
tComplex<T1> operator - (const tComplex<T1>& z, const tComplex<T2>& w)
{
    return tComplex<T1>(z.real - w.real, z.imag - w.imag);
}

template <class T1, class T2> inline
tComplex<T1> operator - (const tComplex<T1>& z)
{
    return tComplex<T1>(-z.real, -z.imag);
}

template <class T1, class T2> inline
tComplex<T1> operator * (const tComplex<T1>& z, const tComplex<T2>& w)
{
    return tComplex<T1>(
        z.real * w.real - z.imag * w.imag,
        z.real * w.imag + z.imag * w.real);
}

template <class T1, class T2> inline
tComplex<T1> operator * (const tComplex<T1>& z, const T2& x)
{
    return tComplex<T1>(x * z.real, x * z.imag);
}

template <class T1, class T2> inline
tComplex<T1> operator * (const T2& x, const tComplex<T1>& z)
{
    return tComplex<T1>(x * z.real, x * z.imag);
}

template <class T1, class T2> inline
tComplex<T1> operator / (const tComplex<T1>& z, const tComplex<T2>& w)
{
    const T1 d = w.real * w.real + w.imag * w.imag;

    return tComplex<T1>(
        (z.real * w.real + z.imag * w.imag) / d,
        (z.imag * w.real - z.real * w.imag) / d);
}

template <class T1, class T2> inline
tComplex<T1> operator / (const tComplex<T1>& z, const T2& x)
{
    return tComplex<T1>(z.real / x, z.imag / x);
}

template <class T> inline
std::ostream& operator << (std::ostream& os, const tComplex<T>& z)
{
  os << "[" << z.real << ", " << z.imag << "]";
  return os;
}

typedef tComplex<float> fComplex;
typedef tComplex<double> dComplex;

#endif
