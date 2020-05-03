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
#ifndef SPHERICAL_HARMONICS_H_
#define SPHERICAL_HARMONICS_H_

#include <cstddef>

#define PVT(l,m) ((m)+((l)*((l)+1))/2)
#define YVR(l,m) ((m)+(l)+((l)*(l)))

namespace au {
namespace edu {
namespace anu {
namespace qm {
namespace ro {
class SphericalHarmonics {
public:
	explicit SphericalHarmonics(int LL);
	~SphericalHarmonics();

	/* Compute an entire set of Y_{l,m}(\theta,\phi) and store in array Y */
	void computeY(const size_t L, const double X, const double phi,
			double* const Y);

private:
	double* A;  // coefficients A for P_l^m
	double* B;  // coefficients B for P_l^m
	double* P;

	/*
	 * Precompute coefficients a_l^m and b_l^m for all l <= L, m <= l
	 */
	void initAB(const size_t LL);

	/* Compute an entire set of P_l^m(x) and store in the array P */
	void computeP(const size_t L, const double X);
};
}  // namespace ro
}  // namespace qm
}  // namespace anu
}  // namespace edu
}  // namespace au

#endif  // SPHERICAL_HARMONICS_H_
