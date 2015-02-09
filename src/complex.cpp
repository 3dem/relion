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
#include "src/complex.h"

// Constructor with two arguments
Complex::Complex(double _r, double _i)
{
    real = _r;
    imag = _i;
}

Complex Complex::operator+ (Complex &op)
{
    return Complex(real + op.real, imag + op.imag);
}

void Complex::operator+= (Complex &op)
{
	real += op.real;
	imag += op.imag;
}

Complex Complex::operator- (Complex &op)
{
    return Complex(real - op.real, imag - op.imag);
}
void Complex::operator-= (Complex &op)
{
	real -= op.real;
	imag -= op.imag;
}

Complex Complex::operator* (Complex &op)
{
    return Complex((real * op.real) - (imag * op.imag), (real * op.imag) + (imag * op.real));
}

Complex Complex::operator* (double op)
{
    return Complex(real*op, imag*op);
}

void Complex::operator*= (double op)
{
	real *= op;
	imag *= op;
}

Complex Complex::operator/(double op)
{
    return Complex(real/op, imag/op);
}

Complex Complex::operator/(Complex &op)
{
    double cd = op.norm();
    double realval = real*op.real + imag*op.imag;
    double imagval = imag*op.real - real*op.imag;
	return Complex(realval/cd, imagval/cd);
}

void Complex::operator/=(double op)
{
	real /= op;
	imag /= op;
}


Complex operator+(const Complex& lhs, const Complex& rhs)
{
	return Complex(lhs.real + rhs.real, lhs.imag + rhs.imag);
}

Complex operator-(const Complex& lhs, const Complex& rhs)
{
	return Complex(lhs.real - rhs.real, lhs.imag - rhs.imag);

}

Complex operator*(const Complex& lhs, const Complex& rhs)
{
	return Complex((lhs.real * rhs.real) - (lhs.imag * rhs.imag), (lhs.real * rhs.imag) + (lhs.imag * rhs.real));
}

Complex operator*(const Complex& lhs, const double& val)
{
	return Complex(lhs.real * val , lhs.imag * val);
}

Complex operator*(const double& val, const Complex& rhs)
{
	return Complex(rhs.real * val , rhs.imag * val);
}

Complex operator/(const Complex& lhs, const double& val)
{
	return Complex(lhs.real / val , lhs.imag / val);
}

void operator+=(Complex& lhs, const Complex& rhs)
{
	lhs.real += rhs.real;
	lhs.imag += rhs.imag;
}
void operator-=(Complex& lhs, const Complex& rhs)
{
	lhs.real -= rhs.real;
	lhs.imag -= rhs.imag;
}

Complex Complex::conj()
{
    return Complex(real, -imag);
}
Complex conj(const Complex& op)
{
	return Complex(op.real, -op.imag);
}


double Complex::abs()
{
    return sqrt(real*real + imag*imag);
}
double abs(const Complex& op)
{
	return sqrt(op.real*op.real + op.imag*op.imag);
}

double Complex::norm()
{
    return real*real + imag*imag;
}
double norm(const Complex& op)
{
	return op.real*op.real + op.imag*op.imag;
}

double Complex::arg()
{
    return atan2(imag, real);
}

double arg(const Complex& op)
{
	return atan2(op.imag, op.real);
}
