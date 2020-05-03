/*  
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright IBM Corporation 2014.
 * (C) Copyright Mahidol University International College 2014.
 */
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "SphericalHarmonics.h"

using au::edu::anu::qm::ro::SphericalHarmonics;

SphericalHarmonics::SphericalHarmonics(int LL) {
	size_t sizeP = (LL + 1) * (LL + 2) / 2;
	P = new double[sizeP];

	A = new double[sizeP];
	B = new double[sizeP];
	this->initAB(LL);
}

SphericalHarmonics::~SphericalHarmonics() {
	delete[] P;
	delete[] A;
	delete[] B;
}

/* 
 * Precompute coefficients a_l^m and b_l^m for all l <= L, m <= l
 */
void SphericalHarmonics::initAB(const size_t LL) {
	assert(A != NULL);
	assert(B != NULL);
	for (size_t l = 2; l <= LL; l++) {
		double ls = l * l;
		double lm1s = (l - 1) * (l - 1);
		for (size_t m = 0; m < l - 1; m++) {
			double ms = m * m;
			A[PVT(l, m)] = sqrt((4 * ls - 1.) / (ls - ms));
			B[PVT(l, m)] = -sqrt((lm1s - ms) / (4 * lm1s - 1.));
		}
	}
}

/* Compute an entire set of P_l^m(x) and store in the array P */
void SphericalHarmonics::computeP(const size_t L, const double X) {
	const double sintheta = sqrt(1. - X * X);
	double temp = 0.39894228040143267794;  // = sqrt(0.5/M_PI)
	P[PVT(0, 0)] = temp;

	if (L > 0) {
		const double SQRT3 = 1.7320508075688772935;
		P[PVT(1, 0)] = X * SQRT3 * temp;
		const double SQRT3DIV2 = -1.2247448713915890491;
		temp = SQRT3DIV2 * sintheta * temp;
		P[PVT(1, 1)] = temp;

		for (size_t l = 2; l <= L; l++) {
			for (size_t m = 0; m < l - 1; m++) {
				P[PVT(l, m)] = A[PVT(l, m)]
						* (X * P[PVT(l - 1, m)]
								+ B[PVT(l, m)] * P[PVT(l - 2, m)]);
			}
			P[PVT(l, l - 1)] = X * sqrt(2 * (l - 1) + 3) * temp;
			temp = -sqrt(1.0 + 0.5 / l) * sintheta * temp;
			P[PVT(l, l)] = temp;
		}
	}
}

/* Compute an entire set of Y_{l,m}(\theta,\phi) and store in the array Y */
void SphericalHarmonics::computeY(const size_t L, const double X,
		const double phi, double* const Y) {
	computeP(L, X);

	for (size_t l = 0; l <= L; l++)
		Y[YVR(l, 0)] = P[PVT(l, 0)] * 0.5 * M_SQRT2;

	// NR2 5.5.4-5.5.5
	double c1 = 1.0, c2 = cos(phi);
	double s1 = 0.0, s2 = -sin(phi);
	double tc = 2.0 * c2;
	for (size_t m = 1; m <= L; m++) {
		double s = tc * s1 - s2;
		double c = tc * c1 - c2;
		s2 = s1;
		s1 = s;
		c2 = c1;
		c1 = c;
		for (size_t l = m; l <= L; l++) {
			Y[YVR(l, -m)] = P[PVT(l, m)] * s;
			Y[YVR(l, m)] = P[PVT(l, m)] * c;
		}
	}
}
