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
#include "src/macros.h"

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


        tComplex& operator += (const tComplex& arg)
        {
            real += arg.real;
            imag += arg.imag;

            return *this;
        }

        tComplex& operator -= (const tComplex& arg)
        {
            real -= arg.real;
            imag -= arg.imag;

            return *this;
        }

        tComplex& operator *= (const tComplex& arg)
        {
            T re = real*arg.real - imag*arg.imag;
            T im = real*arg.imag + imag*arg.real;

            real = re;
            imag = im;

            return *this;
        }

        tComplex& operator /= (const tComplex& arg)
        {
            T cd = arg.real*arg.real + arg.imag*arg.imag;

            T re = real*arg.real + imag*arg.imag;
            T im = imag*arg.real - real*arg.imag;

            real = re/cd;
            imag = im/cd;

            return *this;
        }

        tComplex conj() const
        {
            return tComplex(real, -imag);
        }

        T abs() const
        {
            return sqrt(real*real + imag*imag);
        }

        T norm() const
        {
            return real*real + imag*imag;
        }

        T arg() const
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

template <class T> inline
tComplex<T> operator + (const tComplex<T>& z, const tComplex<T>& w)
{
    return tComplex<T>(z.real + w.real, z.imag + w.imag);
}

template <class T> inline
tComplex<T> operator - (const tComplex<T>& z, const tComplex<T>& w)
{
    return tComplex<T>(z.real - w.real, z.imag - w.imag);
}

template <class T> inline
tComplex<T> operator * (const tComplex<T>& z, const tComplex<T>& w)
{
    return tComplex<T>(
        z.real * w.real - z.imag * w.imag,
        z.real * w.imag + z.imag * w.real);
}

template <class T> inline
tComplex<T> operator * (const tComplex<T>& z, const T& x)
{
    return tComplex<T>(x * z.real, x * z.imag);
}

template <class T> inline
tComplex<T> operator * (const T& x, const tComplex<T>& z)
{
    return tComplex<T>(x * z.real, x * z.imag);
}

template <class T> inline
tComplex<T> operator / (const tComplex<T>& z, const tComplex<T>& w)
{
    const T d = w.real * w.real + w.imag * w.imag;

    return tComplex<T>(
        (z.real * w.real + z.imag * w.imag) / d,
        (z.imag * w.real - z.real * w.imag) / d);
}

template <class T> inline
tComplex<T> operator / (const tComplex<T>& z, const T& x)
{
    return tComplex<T>(z.real / x, z.imag / x);
}

typedef tComplex<float> fComplex;
typedef tComplex<double> dComplex;

#endif
