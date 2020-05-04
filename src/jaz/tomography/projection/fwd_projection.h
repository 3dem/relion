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

#ifndef FORWARD_PROJECTION_H
#define FORWARD_PROJECTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <iostream>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/tomo_stack.h>


class ForwardProjection
{
	public:
		
		template <typename T>
		static void forwardProject(
				const RawImage<tComplex<T>>& reference,	
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<tComplex<T>>& destData,
				RawImage<tComplex<T>>& destPsf,
				double paddingFactor,
				int num_threads = 1);
};
						
						
template <typename T>
void ForwardProjection::forwardProject(
		const RawImage<tComplex<T>>& reference,	
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<T>>& destData,
		RawImage<tComplex<T>>& destPsf,
		double paddingFactor,
		int num_threads)
{
	const int wh2 = destData.xdim;
	const int h2 = destData.ydim;
	const int fc = destData.zdim;
	
	std::vector<gravis::d3Matrix> projTransp(fc);
	
	for (int f = 0; f < fc; f++)
    {
		gravis::d3Matrix A(proj[f](0,0), proj[f](1,0), proj[f](2,0), 
						   proj[f](0,1), proj[f](1,1), proj[f](2,1), 
						   proj[f](0,2), proj[f](1,2), proj[f](2,2) );
				
		projTransp[f] = (1.0 / paddingFactor) * A;
	}
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (long int yi = 0; yi < h2;  yi++)
		for (long int xi = 0; xi < wh2; xi++)
		{
			const double xd = xi;
			const double yd = yi > h2/2? yi - h2 : yi;
					
			const gravis::d3Vector pi(xd, yd, 0.0);
			gravis::d3Vector pw = projTransp[f] * pi;
			
			destData(xi, yi, f) = Interpolation::linearXYZ_FftwHalf_complex(
						reference, pw.x, pw.y, pw.z);
			
			const double ax = std::abs(pw.x);
			const double ay = std::abs(pw.y);
			const double az = std::abs(pw.z);
			
			if (ax < 1.0 && ay < 1.0 && az < 1.0)
			{
				destPsf(xi, yi, f) = (1.0 - ax) * (1.0 - ay) * (1.0 - az);
			}
			else
			{
				destPsf(xi, yi, f) = 0.0;
			}
		}
	}
}

#endif
