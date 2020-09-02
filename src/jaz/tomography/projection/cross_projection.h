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

#ifndef CROSS_PROJECTION_H
#define CROSS_PROJECTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <iostream>

#include "cross_projection.h"

class CrossProjection
{	
	public:
		
		template <typename SrcType, typename DestType>
		static void crossproject_bwd(
			const RawImage<tComplex<SrcType>>& stackFS,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destPSF,
			RawImage<DestType>& destCTF,
			double paddingFactor,
			int num_threads = 1);
};

template <typename SrcType, typename DestType>
void CrossProjection::crossproject_bwd(
		const RawImage<tComplex<SrcType>>& stackFS,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<DestType>>& destFS,
		RawImage<DestType>& destPSF,
		RawImage<DestType>& destCTF,
		double paddingFactor,
		int num_threads)
{
	const int wh = stackFS.xdim;
	const int w = (wh - 1)*2;
	const int h = stackFS.ydim;
	const int hh = h/2 + 1;
	const int fc = stackFS.zdim;
		
	std::vector<gravis::d3Matrix> cubeToSquare(fc), squareToCube(fc);
	
	for (int f = 0; f < fc; f++)
	{
		gravis::d3Matrix A(proj[f](0,0), proj[f](0,1), proj[f](0,2), 
						   proj[f](1,0), proj[f](1,1), proj[f](1,2), 
						   proj[f](2,0), proj[f](2,1), proj[f](2,2) );
				
		cubeToSquare[f] = paddingFactor * A.invert().transpose();
		squareToCube[f] = cubeToSquare[f];
		squareToCube[f].invert();
	}
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int f_dest = 0; f_dest < fc; f_dest++)
	for (long int yi = 0; yi < h; yi++)
	for (long int xi = 0; xi < wh; xi++)
	{
		const double x = xi;
		const double y = yi < h/2? yi : yi - h;
				
		tComplex<DestType> sum(0.0, 0.0);
		DestType wgh = 0.0;
		DestType psf = 0.0;
	
		const gravis::d3Vector p_dest(x,y,0.0);
		const gravis::d3Vector p_3D = squareToCube[f_dest] * p_dest;
	
		for (int f_src = 0; f_src < fc; f_src++)
		{
			if (f_src == f_dest) continue;
			
			const gravis::d3Vector p_src = cubeToSquare[f_src] * p_3D;
			
			if (p_src.z > -1.0 && p_src.z < 1.0 &&
			    std::abs(p_src.x) < wh && std::abs(p_src.y) < hh)
			{
				tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf_clip(stackFS, p_src.x, p_src.y, f_src);
				
				const double cx = 1.0 - std::abs(p_src.x);
				const double cy = 1.0 - std::abs(p_src.y);
				const double cz = 1.0 - std::abs(p_src.z);
				
				sum.real += cz * z0.real;
				sum.imag += cz * z0.imag;
				wgh += cz;
				
				if (cx > 0.0 && cy > 0.0)
				{
					psf += cx * cy * cz;
				}
			}
		}
	
		if (wgh > 0.0) 
		{
			destFS( xi, yi, f_dest) += sum;
			destCTF(xi, yi, f_dest) += wgh;
			destPSF(xi, yi, f_dest) += psf;
		}
	}
}

#endif

