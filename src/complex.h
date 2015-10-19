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
#ifndef COMPLEX_H_
#define COMPLEX_H_
#include <iostream>
#include <cmath>
#include "src/macros.h"

class Complex
{

	public:

	RFLOAT real;
	RFLOAT imag;

    // Constructor
	Complex(RFLOAT _r = 0.0, RFLOAT _i = 0.0);

    Complex operator+(Complex &op);
    void operator+=(Complex &op);

    Complex operator-(Complex &op);
    void operator-=(Complex &op);

    Complex operator*(Complex &op);

    void operator*=(RFLOAT op);

    Complex operator*(RFLOAT op);

    Complex operator/(Complex &op);

    Complex operator/(RFLOAT op);

    void operator/=(RFLOAT op);

    // Complex conjugated
    Complex conj();

    // Abs value: sqrt(real*real+imag*imag)
    RFLOAT abs();

    // Norm value: real*real+imag*imag
    RFLOAT norm();

    // Phase angle: atan2(imag,real)
    RFLOAT arg();


};

Complex conj(const Complex& op);
RFLOAT abs(const Complex& op);
RFLOAT norm(const Complex& op);
RFLOAT arg(const Complex& op);

Complex operator+(const Complex& lhs, const Complex& rhs);
Complex operator-(const Complex& lhs, const Complex& rhs);
Complex operator*(const Complex& lhs, const Complex& rhs);
Complex operator*(const Complex& lhs, const RFLOAT& val);
Complex operator*(const RFLOAT& val, const Complex& rhs);
Complex operator/(const Complex& lhs, const RFLOAT& val);

void operator+=(Complex& lhs, const Complex& rhs);
void operator-=(Complex& lhs, const Complex& rhs);

#endif
