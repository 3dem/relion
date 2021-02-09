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

#include "nelder_mead.h"
#include <src/jaz/util/index_sort.h>
#include <iostream>
#include <cmath>

std::vector<double> NelderMead::optimize(
		const std::vector<double>& initial,
		const Optimization& opt,
		double initialStep, double tolerance, long maxIters,
		double alpha, double gamma, double rho, double sigma,
		bool verbose, double* minCost)
{
	const double n = initial.size();
	const double m = initial.size() + 1;

	std::vector<std::vector<double> > simplex(m);

	simplex[0] = initial;

	for (int j = 1; j < m; j++)
	{
		simplex[j] = initial;
		simplex[j][j-1] += initialStep;
	}

	// avoid runtime allocations (std::vectors live on the heap)
	std::vector<std::vector<double> > nextSimplex(m);
	std::vector<double> values(m), nextValues(m), centroid(n),
			reflected(n), expanded(n), contracted(n);

	void* tempStorage = opt.allocateTempStorage();

	// compute values
	for (int j = 0; j < m; j++)
	{
		values[j] = opt.f(simplex[j], tempStorage);
	}

	for (long i = 0; i < maxIters; i++)
	{
		// sort x and f(x) by ascending f(x)
		std::vector<int> order = IndexSort<double>::sortIndices(values);

		if (verbose)
		{
			opt.report(i, values[order[0]], simplex[order[0]]);
		}

		for (int j = 0; j < m; j++)
		{
			nextSimplex[j] = simplex[order[j]];
			nextValues[j] = values[order[j]];
		}

		simplex = nextSimplex;
		values = nextValues;

		// compute centroid
		for (int k = 0; k < n; k++)
		{
			centroid[k] = 0.0;
		}
		for (int j = 0; j < n; j++) // leave out the worst x
		{
			for (int k = 0; k < n; k++)
			{
				centroid[k] += simplex[j][k];
			}
		}
		for (int k = 0; k < n; k++)
		{
			centroid[k] /= n;
		}

		// check for convergence
		bool allInside = true;
		for (int j = 0; j < m; j++)
		{
			double dx = 0.0;

			for (int k = 0; k < n; k++)
			{
				double ddx = simplex[j][k] - centroid[k];
				dx += ddx * ddx;
			}

			if (sqrt(dx) > tolerance)
			{
				allInside = false;
				break;
			}
		}
		if (allInside)
		{
			if (verbose) std::cout << std::endl;
			return simplex[0];
		}

		// reflect
		for (int k = 0; k < n; k++)
		{
			reflected[k] = (1.0 + alpha) * centroid[k] - alpha * simplex[n][k];
		}

		double vRefl = opt.f(reflected, tempStorage);

		if (vRefl < values[n-1] && vRefl > values[0])
		{
			simplex[n] = reflected;
			values[n] = vRefl;
			continue;
		}

		// expand
		if (vRefl < values[0])
		{
			for (int k = 0; k < n; k++)
			{
				expanded[k] = (1.0 - gamma) * centroid[k] + gamma * reflected[k];
			}

			double vExp = opt.f(expanded, tempStorage);

			if (vExp < vRefl)
			{
				simplex[n] = expanded;
				values[n] = vExp;
			}
			else
			{
				simplex[n] = reflected;
				values[n] = vRefl;
			}

			continue;
		}

		// contract
		for (int k = 0; k < n; k++)
		{
			contracted[k] = (1.0 - rho) * centroid[k] + rho * simplex[n][k];
		}

		double vContr = opt.f(contracted, tempStorage);

		if (vContr < values[n])
		{
			simplex[n] = contracted;
			values[n] = vContr;

			continue;
		}

		// shrink
		for (int j = 1; j < m; j++)
		{
			for (int k = 0; k < n; k++)
			{
				simplex[j][k] = (1.0 - sigma) * simplex[0][k] + sigma * simplex[j][k];
			}

			values[j] = opt.f(simplex[j], tempStorage);
		}
	}

	if (verbose) std::cout << std::endl;

	opt.deallocateTempStorage(tempStorage);

	std::vector<int> order = IndexSort<double>::sortIndices(values);

	if (minCost)
	{
		*minCost = values[order[0]];
	}

	return simplex[order[0]];
}

void NelderMead::test()
{
	RosenbrockBanana rb;
	std::vector<double> initial(2);
	initial[0] = 3.0;
	initial[0] = 1.0;

	std::vector<double> x0 = optimize(initial, rb, 0.5, 0.001, 1000);

	std::cout << "should be close to 1, 1: " << x0[0] << ", " << x0[1] << "\n";
}
