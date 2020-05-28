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

#ifndef FOURIER_BACKPROJECTION_H
#define FOURIER_BACKPROJECTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <iostream>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_voxel.h>


class FourierBackprojection
{
	public:
		
		template <typename SrcType, typename DestType>
		static void backproject_bwd(
			RawImage<tComplex<SrcType>>& stackFS,
			RawImage<SrcType>& weightStack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destPSF,
			RawImage<DestType>& destCTF,
			int num_threads = 1);
		
		template <typename SrcType, typename DestType>
		static void backproject_bwd(
			RawImage<tComplex<SrcType>>& stackFS,
			RawImage<SrcType>& weightStack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destPSF,
			RawImage<DestType>& destCTF,
			RawImage<DestType>& destMP,
			int num_threads = 1);
		
		
		template <typename SrcType, typename DestType>
		static void backprojectSlice_noSF(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			int num_threads);	
		
		template <typename SrcType, typename DestType>
		static void backprojectSlice_noSF(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			RawImage<DestType>& destMP,
			int num_threads);		

		template <typename SrcType, typename DestType>
		static void backprojectSlice_noSF_dualContrast(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& sin_gamma,
			const RawImage<SrcType>& cos_gamma,
			const gravis::d4Matrix& proj,
			RawImage<DualContrastVoxel<DestType>>& dest,
			int num_threads);
		
		
		template <typename DestType>
		static void backprojectSpreadingFunction(
			const gravis::d4Matrix& proj,
			RawImage<DestType>& destPSF);
};


template <typename SrcType, typename DestType>
void FourierBackprojection::backproject_bwd(
		RawImage<tComplex<SrcType>>& stackFS,
		RawImage<SrcType>& weightStack,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<DestType>>& destFS,
		RawImage<DestType>& destPSF,
		RawImage<DestType>& destCTF,
		int num_threads)
{
	const int fc = stackFS.zdim;
	
	for (int f = 0; f < fc; f++)
	{	
		backprojectSlice_noSF(
			stackFS.getSliceRef(f),
			weightStack.getSliceRef(f),
			proj[f],
			destFS,
			destCTF,
			num_threads);
				
		backprojectSpreadingFunction(
			proj[f],
			destPSF);
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backproject_bwd(
		RawImage<tComplex<SrcType>>& stackFS,
		RawImage<SrcType>& weightStack,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<DestType>>& destFS,
		RawImage<DestType>& destPSF,
		RawImage<DestType>& destCTF,
		RawImage<DestType>& destMP,
		int num_threads)
{
	const int fc = stackFS.zdim;
	
	for (int f = 0; f < fc; f++)
	{	
		backprojectSlice_noSF(
			stackFS.getSliceRef(f),
			weightStack.getSliceRef(f),
			proj[f],
			destFS,
			destCTF,
			destMP,
			num_threads);
				
		backprojectSpreadingFunction(proj[f], destPSF);
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_noSF(
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF,
				int num_threads)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;
	
	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_noSF: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	gravis::d3Matrix A(proj(0,0), proj(0,1), proj(0,2), 
					   proj(1,0), proj(1,1), proj(1,2), 
					   proj(2,0), proj(2,1), proj(2,2) );
			
	gravis::d3Matrix projInvTransp = A.invert().transpose();
	gravis::d3Vector normal(projInvTransp(2,0), projInvTransp(2,1), projInvTransp(2,2));
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int z = 0; z < d3; z++)
	for (long int y = 0; y < h3; y++)
	{
		const double yy = y > h3/2? y - h3 : y;
		const double zz = z > d3/2? z - d3 : z;
		
		const double yz = normal.y * yy + normal.z * zz;
		
		long int x0, x1;
		
		if (normal.x == 0.0) 
		{
			if (yz > -1.0 && yz < 1.0)
			{
				x0 = 0;
				x1 = wh3-1;
			}
			else
			{
				x0 = 0;
				x1 = -1;
			}
		}
		else
		{
			const double a0 = (-yz - 1.0) / normal.x;
			const double a1 = (-yz + 1.0) / normal.x;
			
			if (a0 < a1)
			{
				x0 = std::ceil(a0);
				x1 = std::floor(a1);
			}
			else
			{
				x0 = std::ceil(a1);
				x1 = std::floor(a0);
			}
			
			if (x0 < 0) x0 = 0;
			if (x1 > wh3-1) x1 = wh3-1;
		}
		
		for (long int x = x0; x <= x1; x++)
		{
			gravis::d3Vector pw(x,yy,zz);		
			gravis::d3Vector pi = projInvTransp * pw;
			
			if (pi.z > -1.0 && pi.z < 1.0 &&
				std::abs(pi.x) < wh2 && std::abs(pi.y) < h2/2 + 1 )
			{
				const double c = 1.0 - std::abs(pi.z);
				
				tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf(dataFS, pi.x, pi.y, 0);
				const DestType wgh = Interpolation::linearXY_symmetric_FftwHalf(weight, pi.x, pi.y, 0);
								
				destFS(x,y,z) += tComplex<DestType>(c * z0.real, c * z0.imag);
				destCTF(x,y,z) += c * wgh;
			}
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_noSF(
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF,
				RawImage<DestType>& destMP,
				int num_threads)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;
	
	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_noSF: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	gravis::d3Matrix A(proj(0,0), proj(0,1), proj(0,2), 
					   proj(1,0), proj(1,1), proj(1,2), 
					   proj(2,0), proj(2,1), proj(2,2) );
			
	gravis::d3Matrix projInvTransp = A.invert().transpose();
	gravis::d3Vector normal(projInvTransp(2,0), projInvTransp(2,1), projInvTransp(2,2));
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int z = 0; z < d3; z++)
	for (long int y = 0; y < h3; y++)
	{
		const double yy = y > h3/2? y - h3 : y;
		const double zz = z > d3/2? z - d3 : z;
		
		const double yz = normal.y * yy + normal.z * zz;
		
		long int x0, x1;
		
		if (normal.x == 0.0) 
		{
			if (yz > -1.0 && yz < 1.0)
			{
				x0 = 0;
				x1 = wh3-1;
			}
			else
			{
				x0 = 0;
				x1 = -1;
			}
		}
		else
		{
			const double a0 = (-yz - 1.0) / normal.x;
			const double a1 = (-yz + 1.0) / normal.x;
			
			if (a0 < a1)
			{
				x0 = std::ceil(a0);
				x1 = std::floor(a1);
			}
			else
			{
				x0 = std::ceil(a1);
				x1 = std::floor(a0);
			}
			
			if (x0 < 0) x0 = 0;
			if (x1 > wh3-1) x1 = wh3-1;
		}
		
		for (long int x = x0; x <= x1; x++)
		{
			gravis::d3Vector pw(x,yy,zz);		
			gravis::d3Vector pi = projInvTransp * pw;
			
			if (pi.z > -1.0 && pi.z < 1.0 &&
				std::abs(pi.x) < wh2 && std::abs(pi.y) < h2/2 + 1 )
			{
				const double c = 1.0 - std::abs(pi.z);
				
				tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf(dataFS, pi.x, pi.y, 0);
				const DestType wgh = Interpolation::linearXY_symmetric_FftwHalf(weight, pi.x, pi.y, 0);
								
				destFS(x,y,z) += tComplex<DestType>(c * z0.real, c * z0.imag);
				destCTF(x,y,z) += c * wgh;
				destMP(x,y,z) += c;
			}
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_noSF_dualContrast(
	const RawImage<tComplex<SrcType>>& dataFS,
	const RawImage<SrcType>& sin_gamma,
	const RawImage<SrcType>& cos_gamma,
	const gravis::d4Matrix& proj,
	RawImage<DualContrastVoxel<DestType>>& dest,
	int num_threads)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;

	const int wh3 = dest.xdim;
	const int h3 = dest.ydim;
	const int d3 = dest.zdim;

	gravis::d3Matrix A(proj(0,0), proj(0,1), proj(0,2),
					   proj(1,0), proj(1,1), proj(1,2),
					   proj(2,0), proj(2,1), proj(2,2) );

	gravis::d3Matrix projInvTransp = A.invert().transpose();
	gravis::d3Vector normal(projInvTransp(2,0), projInvTransp(2,1), projInvTransp(2,2));

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d3; z++)
	for (long int y = 0; y < h3; y++)
	{
		const double yy = y > h3/2? y - h3 : y;
		const double zz = z > d3/2? z - d3 : z;

		const double yz = normal.y * yy + normal.z * zz;

		long int x0, x1;

		if (normal.x == 0.0)
		{
			if (yz > -1.0 && yz < 1.0)
			{
				x0 = 0;
				x1 = wh3-1;
			}
			else
			{
				x0 = 0;
				x1 = -1;
			}
		}
		else
		{
			const double a0 = (-yz - 1.0) / normal.x;
			const double a1 = (-yz + 1.0) / normal.x;

			if (a0 < a1)
			{
				x0 = std::ceil(a0);
				x1 = std::floor(a1);
			}
			else
			{
				x0 = std::ceil(a1);
				x1 = std::floor(a0);
			}

			if (x0 < 0) x0 = 0;
			if (x1 > wh3-1) x1 = wh3-1;
		}

		for (long int x = x0; x <= x1; x++)
		{
			gravis::d3Vector pw(x,yy,zz);
			gravis::d3Vector pi = projInvTransp * pw;

			if (pi.z > -1.0 && pi.z < 1.0 &&
				std::abs(pi.x) < wh2 && std::abs(pi.y) < h2/2 + 1 )
			{
				const double c = 1.0 - std::abs(pi.z);

				const tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf(dataFS, pi.x, pi.y, 0);
				const tComplex<DestType> z1(z0.real, z0.imag);

				const DestType sin_g = Interpolation::linearXY_symmetric_FftwHalf(sin_gamma, pi.x, pi.y, 0);
				const DestType cos_g = Interpolation::linearXY_symmetric_FftwHalf(cos_gamma, pi.x, pi.y, 0);

				dest(x,y,z).data_sin += c * sin_g * z1;
				dest(x,y,z).data_cos += c * cos_g * z1;

				dest(x,y,z).weight_sin2    += c * sin_g * sin_g;
				dest(x,y,z).weight_sin_cos += c * sin_g * cos_g;
				dest(x,y,z).weight_cos2    += c * cos_g * cos_g;
			}
		}
	}
}

template <typename DestType>
void FourierBackprojection::backprojectSpreadingFunction(
	const gravis::d4Matrix& proj,
	RawImage<DestType>& destPSF)
{	
	const int h3 = destPSF.ydim;
	const int d3 = destPSF.zdim;
	
	gravis::d3Matrix A(proj(0,0), proj(0,1), proj(0,2), 
					   proj(1,0), proj(1,1), proj(1,2), 
					   proj(2,0), proj(2,1), proj(2,2) );
			
	gravis::d3Matrix projInvTransp = A.invert().transpose();
	
	const double rng = std::floor( (1.0 / sqrt(std::abs(projInvTransp.det())) ) + 1e-4);
	
	for (int dz = -rng; dz <= rng; dz++)
	for (int dy = -rng; dy <= rng; dy++)
	for (int dx = 0; dx <= rng; dx++)
	{
		gravis::d3Vector pw(dx,dy,dz);		
		gravis::d3Vector pi = projInvTransp * pw;
		
		const double wx = 1.0 - std::abs(pi.x);
		const double wy = 1.0 - std::abs(pi.y);
		const double wz = 1.0 - std::abs(pi.z);
		
		if (wx > 0.0 && wy > 0.0 && wz > 0.0)
		{
			const int x = dx;
			const int y = (h3 + dy) % h3;
			const int z = (d3 + dz) % d3;
			
			destPSF(x,y,z) += wx * wy * wz;
		}		
	}
}

#endif
