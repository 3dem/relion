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

class Complex
{

	public:

	double real;
	double imag;

    // Constructor
	Complex(double _r = 0.0, double _i = 0.0);

    Complex operator+(Complex &op);
    void operator+=(Complex &op);

    Complex operator-(Complex &op);
    void operator-=(Complex &op);

    Complex operator*(Complex &op);

    void operator*=(double op);

    Complex operator*(double op);

    Complex operator/(Complex &op);

    Complex operator/(double op);

    void operator/=(double op);

    // Complex conjugated
    Complex conj();

    // Abs value: sqrt(real*real+imag*imag)
    double abs();

    // Norm value: real*real+imag*imag
    double norm();

    // Phase angle: atan2(imag,real)
    double arg();


};

Complex conj(const Complex& op);
double abs(const Complex& op);
double norm(const Complex& op);
double arg(const Complex& op);

Complex operator+(const Complex& lhs, const Complex& rhs);
Complex operator-(const Complex& lhs, const Complex& rhs);
Complex operator*(const Complex& lhs, const Complex& rhs);
Complex operator*(const Complex& lhs, const double& val);
Complex operator*(const double& val, const Complex& rhs);
Complex operator/(const Complex& lhs, const double& val);

void operator+=(Complex& lhs, const Complex& rhs);
void operator-=(Complex& lhs, const Complex& rhs);

#endif
