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

#include "optimization.h"
#include <iostream>
#include <iomanip>

double RosenbrockBanana::f(const std::vector<double> &x, void *tempStorage) const
{
    const double a = 1.0;
    const double b = 100.0;
    const double& xx = x[0];
    const double& yy = x[1];

    return (a - xx) * (a - xx) + b * (yy - xx * xx) * (yy - xx * xx);
}

void RosenbrockBanana::grad(
    const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
    const double a = 1.0;
    const double b = 100.0;
    const double& xx = x[0];
    const double& yy = x[1];

    gradDest[0] = -2.0 * (a - xx) - 4.0 * b * (yy - xx * xx) * xx;
    gradDest[1] = 2.0 * b * (yy - xx * xx);
}

void DifferentiableOptimization::testGradient(const std::vector<double> &x, double eps)
{
	const int n = x.size();	
	
	std::vector<double> x1 = x, grad0(n,0.0);
	
	const double f0 = f(x, 0);
	grad(x, grad0, 0);
	
	for (int i = 0; i < n; i++)
	{
		x1[i] = x[i] + eps;		
		const double f1 = f(x1, 0);
		
		x1[i] = x[i];
		
		std::cout << i << ": " << std::setprecision(16) << (f1 - f0)/eps << " vs. " << grad0[i] << "\n";		
	}
}
