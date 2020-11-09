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

#ifndef JAZ_OPTIMIZATION_H
#define JAZ_OPTIMIZATION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>


// abstract class for optimization problems
class Optimization
{
	public:
		
		virtual double f(const std::vector<double>& x, void* tempStorage) const = 0;
		
		virtual void* allocateTempStorage() const 
		{
			return 0;
		}
		
		virtual void deallocateTempStorage(void* ts) const 
		{}
		
		virtual void report(int iteration, double cost, const std::vector<double>& x) const 
		{
			std::cout.precision(16);

			if (iteration == 0) 
			{
				return;
			}
			else if (iteration == 1) 
			{
				std::cout << "iteration    cost\n";
				std::cout << ' ' << std::setw(6) << iteration << "       " << cost;
			}
			else
			{ 
				std::cout << "\r                                                                         ";	
				std::cout << '\r' << ' ' << std::setw(6) << iteration << "       " << cost;
			}
			
			std::cout.flush();
		}
};

class DifferentiableOptimization : public Optimization
{
	public:

		virtual void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const = 0;

		void testGradient(const std::vector<double>& x, double eps);
};

class FastDifferentiableOptimization
{
	public:

		virtual double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const = 0;

		void testGradient(const std::vector<double>& x, double eps);

		virtual void report(int iteration, double cost, const std::vector<double>& x) const
		{
			std::cout.precision(16);

			if (iteration == 0)
			{
				return;
			}
			else if (iteration == 1)
			{
				std::cout << "iteration    cost\n";
				std::cout << ' ' << std::setw(6) << iteration << "       " << cost;
			}
			else
			{
				std::cout << "\r                                                                         ";
				std::cout << '\r' << ' ' << std::setw(6) << iteration << "       " << cost;
			}

			std::cout.flush();
		}
};

class RosenbrockBanana : public DifferentiableOptimization
{
	public:
		
		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
};

#endif
