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
		static void backprojectSlice_backward(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			int num_threads);

		template <typename SrcType, typename DestType>
		static void backprojectSlice_backward(
			int maxFreq,
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			int num_threads);

		template <typename SrcType, typename DestType>
		static void backprojectSlice_backward_withMultiplicity(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			RawImage<DestType>& destMP,
			int num_threads);		

		template <typename SrcType, typename DestType>
		static void backprojectSlice_dualContrast_backward(
			const RawImage<tComplex<SrcType>>& sin_gamma_data,
			const RawImage<tComplex<SrcType>>& cos_gamma_data,
			const RawImage<SrcType>& sin2_weight,
			const RawImage<SrcType>& sin_cos_weight,
			const RawImage<SrcType>& cos2_weight,
			const gravis::d4Matrix& proj,
			RawImage<DualContrastVoxel<DestType>>& dest,
			int num_threads);
		
		template <typename SrcType, typename DestType, class PointInsertion>
		static void backprojectSlice_forward(
			const PointInsertion& pointInsertion,
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF);

		template <typename SrcType, typename DestType>
		static void backprojectSlice_forward_with_multiplicity(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			RawImage<DestType>& destMult);

		template <typename SrcType, typename DestType>
		static void backprojectSlice_forward_with_multiplicity(
			const int* xRanges,
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			RawImage<DestType>& destMult);
		
		
		
		template <typename SrcType, typename DestType>
		static void backprojectStack_backward(
			RawImage<tComplex<SrcType>>& stackFS,
			RawImage<SrcType>& weightStack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destPSF,
			RawImage<DestType>& destCTF,
			int num_threads = 1);
		
		template <typename SrcType, typename DestType>
		static void backprojectStack_backward(
			RawImage<tComplex<SrcType>>& stackFS,
			RawImage<SrcType>& weightStack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destPSF,
			RawImage<DestType>& destCTF,
			RawImage<DestType>& destMP,
			int num_threads = 1);
		

		
		template <typename SrcType, typename DestType>
		static inline void rasteriseLine_Ewald(
			const RawImage<tComplex<SrcType>>& dataPQ,
			const RawImage<SrcType>& weight,
			int x0, 
			int x1,
			const int y, 
			const int z,
			const double D,
			const gravis::d3Matrix& A,
			const double yy, 
			const double zz,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF);
		
		
		template <typename SrcType, typename DestType>
		static void backprojectSphere_backward(
			const RawImage<tComplex<SrcType>>& dataPQ,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			const double radius,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			int num_threads);
		
		template <typename SrcType, typename DestType>
		static void backprojectSphere_backward_slow(
			const RawImage<tComplex<SrcType>>& dataPQ,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			const double radius,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			int num_threads);
		
		/*template <typename SrcType, typename DestType>
		static void backprojectSlice_dualContrast_backward_curved_slow(
			const RawImage<tComplex<SrcType>>& dataPQ,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			const double radius,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF,
			int num_threads);*/
		
		

		
		template <typename SrcType, typename DestType, class PointInsertion>
		static void backprojectSphere_forward(
			const PointInsertion& pointInsertion,
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			double radius,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF);
		
		template <typename SrcType, typename DestType, class PointInsertion>
		static void backprojectSphere_dualContrast_forward(
			const PointInsertion& pointInsertion,
			const RawImage<tComplex<SrcType>>& sin_gamma_data,
			const RawImage<tComplex<SrcType>>& cos_gamma_data,
			const RawImage<SrcType>& sin2_weight,
			const RawImage<SrcType>& sin_cos_weight,
			const RawImage<SrcType>& cos2_weight,
			const gravis::d4Matrix& proj,
			double radius,
			bool conjugate,
			RawImage<DualContrastVoxel<DestType>>& dest);
		
		template <typename SrcType, typename DestType>
		static void backprojectSlice_forward_wrap(
			const RawImage<tComplex<SrcType>>& dataFS,
			const RawImage<SrcType>& weight,
			const gravis::d4Matrix& proj,
			RawImage<tComplex<DestType>>& destFS,
			RawImage<DestType>& destCTF);
		
		
		
		template <typename DestType>
		static void backprojectSpreadingFunction(
			const gravis::d4Matrix& proj,
			RawImage<DestType>& destPSF);
};


template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectStack_backward(
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
		backprojectSlice_backward(
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
void FourierBackprojection::backprojectStack_backward(
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
		backprojectSlice_backward_withMultiplicity(
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
void FourierBackprojection::backprojectSlice_backward(
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
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_backward: destCTF has wrong size ("
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
		const double yy = y >= h3/2? y - h3 : y;
		const double zz = z >= d3/2? z - d3 : z;

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

				tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf_clip(dataFS, pi.x, pi.y, 0);
				const DestType wgh = Interpolation::linearXY_symmetric_FftwHalf_clip(weight, pi.x, pi.y, 0);

				destFS(x,y,z) += tComplex<DestType>(c * z0.real, c * z0.imag);
				destCTF(x,y,z) += c * wgh;
			}
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_backward(
				int maxFreq,
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
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_backward: destCTF has wrong size ("
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
		const double yy = y >= h3/2? y - h3 : y;
		const double zz = z >= d3/2? z - d3 : z;

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

		const int max_x = (int) sqrt(maxFreq*maxFreq - yy*yy - zz*zz);

		if (x1 > max_x) x1 = max_x;

		for (long int x = x0; x <= x1; x++)
		{
			gravis::d3Vector pw(x,yy,zz);
			gravis::d3Vector pi = projInvTransp * pw;

			if (pi.z > -1.0 && pi.z < 1.0 &&
				std::abs(pi.x) < wh2 && std::abs(pi.y) < h2/2 + 1 )
			{
				const double c = 1.0 - std::abs(pi.z);

				tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf_clip(dataFS, pi.x, pi.y, 0);
				const DestType wgh = Interpolation::linearXY_symmetric_FftwHalf_clip(weight, pi.x, pi.y, 0);

				destFS(x,y,z) += tComplex<DestType>(c * z0.real, c * z0.imag);
				destCTF(x,y,z) += c * wgh;
			}
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_backward_withMultiplicity(
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
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_backward_withMultiplicity: destCTF has wrong size ("
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
		const double yy = y >= h3/2? y - h3 : y;
		const double zz = z >= d3/2? z - d3 : z;
		
		const double yz = normal.y * yy + normal.z * zz;
		
		const int max_x = wh3-1;
		
		long int x0, x1;
		
		if (normal.x == 0.0) 
		{
			if (yz > -1.0 && yz < 1.0)
			{
				x0 = 0;
				x1 = max_x;
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
			if (x1 > max_x) x1 = max_x;
		}
		
		for (long int x = x0; x <= x1; x++)
		{
			gravis::d3Vector pw(x,yy,zz);		
			gravis::d3Vector pi = projInvTransp * pw;
			
			if (pi.z > -1.0 && pi.z < 1.0 &&
				std::abs(pi.x) < wh2 && std::abs(pi.y) < h2/2 + 1 )
			{
				const double c = 1.0 - std::abs(pi.z);
				
				tComplex<SrcType> z0 = Interpolation::linearXY_complex_FftwHalf_clip(dataFS, pi.x, pi.y, 0);
				const DestType wgh = Interpolation::linearXY_symmetric_FftwHalf_clip(weight, pi.x, pi.y, 0);
								
				destFS(x,y,z) += tComplex<DestType>(c * z0.real, c * z0.imag);
				destCTF(x,y,z) += c * wgh;
				destMP(x,y,z) += c;
			}
		}
	}
}


template <typename SrcType, typename DestType>
inline void FourierBackprojection::rasteriseLine_Ewald(
        const RawImage<tComplex<SrcType>>& dataPQ,
		const RawImage<SrcType>& weight,
        int x0, 
        int x1,
        const int y, 
        const int z,
        const double D,
        const gravis::d3Matrix& A,
        const double yy, 
        const double zz,
        RawImage<tComplex<DestType>>& destFS,
		RawImage<DestType>& destCTF)
{
	const int w2 = dataPQ.xdim;
	const int h2 = dataPQ.ydim;
	const int m = w2 / 2;	
	
	const int wh3 = destFS.xdim;	
	
	if (x0 < 0) x0 = 0;
	if (x1 >= wh3) x1 = wh3 - 1;	
	
	for (long int x = x0; x <= x1; x++)
	{
		const gravis::d3Vector pw(x,yy,zz);		
		const gravis::d3Vector pi0 = A * pw;
		const gravis::d3Vector pi(pi0.x, pi0.y, pi0.z - D * (pi0.x * pi0.x + pi0.y * pi0.y));
		
		const double u = m + pi.x;
		const double v = m + pi.y;
		
		        
		if (pi.z > -1.0 && pi.z < 1.0 &&
			u >= 0 && u <= w2-1 && 
			v >= 0 && v <= h2-1 )
		{
			const double c = 1.0 - std::abs(pi.z);
			
			tComplex<SrcType> z0 = Interpolation::linearXY_clip(dataPQ, u, v, 0);
			const DestType wgh = Interpolation::linearXY_symmetric_FftwHalf_clip(weight, pi.x, pi.y, 0);

			destFS(x,y,z) += tComplex<DestType>(c * z0.real, c * z0.imag);
			destCTF(x,y,z) += c * wgh;
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSphere_backward(
				const RawImage<tComplex<SrcType>>& dataPQ,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				const double radius,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF,
				int num_threads)
{
	const int wh3 = destFS.xdim;
	const int h3  = destFS.ydim;
	const int d3  = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSphere_backward: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	gravis::d3Matrix Pt(proj(0,0), proj(1,0), proj(2,0), 
					   proj(0,1), proj(1,1), proj(2,1), 
					   proj(0,2), proj(1,2), proj(2,2) );
			
	gravis::d3Matrix A = Pt.invert();
	
	gravis::d3Vector ax(A(0,0), A(0,1), A(0,2));
	gravis::d3Vector ay(A(1,0), A(1,1), A(1,2));
	gravis::d3Vector az(A(2,0), A(2,1), A(2,2));
	
	const double D = 1.0 / (2.0 * radius);
		
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int z = 0; z < d3; z++)
	for (long int y = 0; y < h3; y++)
	{
		const double yy = y >= h3/2? y - h3 : y;
		const double zz = z >= d3/2? z - d3 : z;
		
		/*
			forward model:
			
				3x3 Matrix P = [Pu, Pv, Pw]
			
				(u,v) -> (x,y,z) = u * Pu + v * Pv + D * (u² + v²) * Ph
			
			
			extend:
			
				(u,v,h) -> (x,y,z)   
				
				(x,y,z) = u * Pu + v * Pv + [D * (u² + v²) + h] * Ph
						= P (u,v,h)  + D * (u² + v²) * Ph
						=: (x,y,z)_0 + D * (u² + v²) * Ph				

			invert:
			
				Matrix A := P^(-1) = [Ax, Ay, Az] = [Au, Av, Ah]^T
				
				(x,y,z)_0 -> (u,v,h) 
				
				(u,v,h) = A * (x,y,z)_0
				        = A * [(x,y,z) - D * (u² + v²) * Ph]
				        = A * (x,y,z) - D * (u² + v²) * eh
				
				        = [ Au * (x,y,z) ]
				          [ Av * (x,y,z) ]
				          [ Ah * (x,y,z) - D * (u² + v²)]
 
				
				        = [ Au * (x,y,z) ]
				          [ Av * (x,y,z) ]
				          [ Ah * (x,y,z) - D * ((Au * (x,y,z))² + (Av * (x,y,z))²)]

			given y and z, find the range of x for which -1 < h < 1:
			
				Ah * (x,y,z) - D * ((Au * (x,y,z))² + (Av * (x,y,z))²)  in [-1,1]
				
				= Ahx * x + Ahy * y + Ahz * z
				  - D * ((Aux * x + Auy * y + Auz * z)² + (Avx * x + Avy * y + Avz * z)²)
				  
				= Axz * x + Ayz * y + Azz * z
				  - D * ((Axx * x + Ayx * y + Azx * z)² + (Axy * x + Ayy * y + Azy * z)²)
				  
				= Axz * x + Ayz * y + Azz * z				
				  - D * (
				           Axx² x² + 2 Axx Ayx x y + 2 Axx Azx x z 
				         + Ayx² y² + 2 Ayx Azx y z
				         + Azx² z² 
				
				         + Axy² x² + 2 Axy Ayy x y + 2 Axy Azy x z 
				         + Ayy² y² + 2 Ayy Azy y z
				         + Azy² z² )
				

				= -D [Axx² + Axy²] * x²
				
				  + [-D {2(Axx Ayx + Axy Ayy) y + 2(Axx Azx + Axy Azy) z} + Axz] * x
				  
				  + -D {(Ayx²+Ayy²) y² + 2(Ayx Azx + Ayy Azy) y z + (Azx² + Azy²) z²} + Ayz * y + Azz * z
				  
			
				=: f(x) =: alpha x² + beta x + gamma
						  
		*/
		
		
		const double alpha = -D * (ax.x * ax.x + ax.y * ax.y);
		
		const double beta  =
		        -D * (
						 2 * (ax.x * ay.x + ax.y * ay.y) * yy
					   + 2 * (ax.x * az.x + ax.y * az.y) * zz 
		             )
				+ ax.z;
		
		const double gamma_0 = 
		        -D * (      
							 (ay.x * ay.x + ay.y * ay.y) * yy * yy 
					   + 2 * (ay.x * az.x + ay.y * az.y) * yy * zz
					   +     (az.x * az.x + az.y * az.y) * zz * zz 
					 ) 
				+ ay.z * yy + az.z * zz;
		
		const double gamma_up   = alpha > 0.0? gamma_0 + 1 :  gamma_0 - 1;
		const double gamma_down = alpha > 0.0? gamma_0 - 1 :  gamma_0 + 1;
		
		const double discr_out  = beta * beta - 4 * alpha * gamma_down;
		const double discr_in   = beta * beta - 4 * alpha * gamma_up;
		        
		if (discr_out < 0.0) continue;
		
		const double mid = -0.5 * beta / alpha;
		const double outer_step = 0.5 * sqrt(discr_out) / std::abs(alpha);
		
		
		if (discr_in < 0.0)
		{
			const int x0 = std::ceil(mid - outer_step);
			const int x1 = std::floor(mid + outer_step);
			
			rasteriseLine_Ewald(
				dataPQ, weight, x0, x1, y, z, D, A,
				yy, zz, destFS, destCTF);
		}
		else
		{
			const double inner_step = 0.5 * sqrt(discr_in) / std::abs(alpha);
			
			const int x0_0 = std::ceil( mid - outer_step);
			const int x1_0 = std::floor(mid - inner_step);
			
			const int x0_1 = std::ceil( mid + inner_step);
			const int x1_1 = std::floor(mid + outer_step);
						
			rasteriseLine_Ewald(
				dataPQ, weight, x0_0, x1_0, y, z, D, A,
				yy, zz, destFS, destCTF);
			
			rasteriseLine_Ewald(
				dataPQ, weight, x0_1, x1_1, y, z, D, A,
				yy, zz, destFS, destCTF);
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSphere_backward_slow(
				const RawImage<tComplex<SrcType>>& dataPQ,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				const double radius,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF,
				int num_threads)
{
	const int wh3 = destFS.xdim;
	const int h3  = destFS.ydim;
	const int d3  = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSphere_backward_slow: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	gravis::d3Matrix Pt(proj(0,0), proj(1,0), proj(2,0), 
	                    proj(0,1), proj(1,1), proj(2,1), 
	                    proj(0,2), proj(1,2), proj(2,2) );
			
	gravis::d3Matrix A = Pt.invert();
	
	const double D = 1.0 / (2.0 * radius);
		
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int z = 0; z < d3; z++)
	for (long int y = 0; y < h3; y++)
	{
		const double yy = y >= h3/2? y - h3 : y;
		const double zz = z >= d3/2? z - d3 : z;
		
		rasteriseLine_Ewald(
			dataPQ, weight, 0, wh3-1, y, z, D, A,
			yy, zz, destFS, destCTF);
	}
}


template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_dualContrast_backward(
	const RawImage<tComplex<SrcType>>& sin_gamma_data,
	const RawImage<tComplex<SrcType>>& cos_gamma_data,
	const RawImage<SrcType>& sin2_weight,
	const RawImage<SrcType>& sin_cos_weight,
	const RawImage<SrcType>& cos2_weight,
	const gravis::d4Matrix& proj,
	RawImage<DualContrastVoxel<DestType>>& dest,
	int num_threads)
{
	const int wh2 = sin_gamma_data.xdim;
	const int h2 = sin_gamma_data.ydim;

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
		const double yy = y >= h3/2? y - h3 : y;
		const double zz = z >= d3/2? z - d3 : z;

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

				const tComplex<SrcType> zs0 = Interpolation::linearXY_complex_FftwHalf_clip(
							sin_gamma_data, pi.x, pi.y, 0);

				const tComplex<DestType> zs1(zs0.real, zs0.imag);

				const tComplex<SrcType> zc0 = Interpolation::linearXY_complex_FftwHalf_clip(
							cos_gamma_data, pi.x, pi.y, 0);

				const tComplex<DestType> zc1(zc0.real, zc0.imag);

				const DestType sin2_g = Interpolation::linearXY_symmetric_FftwHalf_clip(
							sin2_weight, pi.x, pi.y, 0);

				const DestType sin_cos_g = Interpolation::linearXY_symmetric_FftwHalf_clip(
							sin_cos_weight, pi.x, pi.y, 0);

				const DestType cos2_g = Interpolation::linearXY_symmetric_FftwHalf_clip(
							cos2_weight, pi.x, pi.y, 0);

				dest(x,y,z).data_sin += c * zs1;
				dest(x,y,z).data_cos += c * zc1;

				dest(x,y,z).weight_sin2    += c * sin2_g;
				dest(x,y,z).weight_sin_cos += c * sin_cos_g;
				dest(x,y,z).weight_cos2    += c * cos2_g;
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

template <typename SrcType, typename DestType, class PointInsertion>
void FourierBackprojection::backprojectSlice_forward(
				const PointInsertion& pointInsertion,
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;
	
	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_forward: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	const gravis::d3Vector u(proj(0,0), proj(0,1), proj(0,2));
	const gravis::d3Vector v(proj(1,0), proj(1,1), proj(1,2));
	
	for (long int y = 0; y < h2;  y++)
	for (long int x = (y > 0 && y < h2/2? 0 : 1); x < wh2; x++)
	{
		const double xx = x;
		const double yy = y < h2/2? y : y - h2;
		
		gravis::d3Vector pos3 = xx * u + yy * v;
		
		bool conj = false;
		
		if (pos3.x < 0)
		{
			pos3 = -pos3;
			conj = true;
		}		
		
		const tComplex<SrcType> value = conj? dataFS(x,y).conj() : dataFS(x,y);
		const tComplex<SrcType> wgh = weight(x,y);
		
		pointInsertion.insert(value, wgh, pos3, destFS, destCTF);
	}
}



template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_forward_with_multiplicity(
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF,
				RawImage<DestType>& destMult)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;

	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;

	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_forward_with_multiplicity: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}

	const gravis::d3Vector u(proj(0,0), proj(0,1), proj(0,2));
	const gravis::d3Vector v(proj(1,0), proj(1,1), proj(1,2));
	const double scale = u.length() * v.length();

	for (long int y = 0; y < h2;  y++)
	for (long int x = (y < h2/2? 0 : 1); x < wh2; x++)
	//for (long int x = (y > 0 && y < h2/2? 0 : 1); x < wh2; x++)  // exclude the origin pixel
	{
		const double xx = x;
		const double yy = y < h2/2? y : y - h2;

		gravis::d3Vector pos3 = xx * u + yy * v;

		bool conj = false;

		if (pos3.x < 0)
		{
			pos3 = -pos3;
			conj = true;
		}

		const tComplex<SrcType> value = conj? dataFS(x,y).conj() : dataFS(x,y);
		const tComplex<SrcType> wgh = weight(x,y);

		{
			const int x0 = std::floor(pos3.x);
			const int y0 = std::floor(pos3.y);
			const int z0 = std::floor(pos3.z);

			for (int dz = 0; dz < 2; dz++)
			for (int dy = 0; dy < 2; dy++)
			for (int dx = 0; dx < 2; dx++)
			{
				const int xg = x0 + dx;
				const int yg = y0 + dy;
				const int zg = z0 + dz;

				if ( xg < wh3
				  && yg >= -h3/2 && yg < h3/2
				  && zg >= -d3/2 && zg < d3/2)
				{
					const int xi = xg;
					const int yi = yg >= 0? yg : yg + h3;
					const int zi = zg >= 0? zg : zg + d3;

					const double fx = 1.0 - std::abs(pos3.x - xg);
					const double fy = 1.0 - std::abs(pos3.y - yg);
					const double fz = 1.0 - std::abs(pos3.z - zg);

					const double m = scale * fx * fy * fz;

					destFS(  xi,yi,zi) += m * value;
					destCTF( xi,yi,zi) += m * wgh;
					destMult(xi,yi,zi) += m;

					if (xi == 0)
					{
						const int yim = (h3 - yi) % h3;
						const int zim = (d3 - zi) % d3;

						destFS(  0,yim,zim) += m * value.conj();
						destCTF( 0,yim,zim) += m * wgh;
						destMult(0,yim,zim) += m;
					}
				}
			}
		}
	}
}

template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_forward_with_multiplicity(
				const int* xRanges,
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF,
				RawImage<DestType>& destMult)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;

	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;

	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_forward_with_multiplicity: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}

	const gravis::d3Vector u(proj(0,0), proj(0,1), proj(0,2));
	const gravis::d3Vector v(proj(1,0), proj(1,1), proj(1,2));
	const double scale = u.length() * v.length();

	for (long int y = 0; y < h2;  y++)
	{
		for (long int x = (y < h2/2? 0 : 1); x < xRanges[y]; x++)
		{
			const double xx = x;
			const double yy = y < h2/2? y : y - h2;

			gravis::d3Vector pos3 = xx * u + yy * v;

			bool conj = false;

			if (pos3.x < 0)
			{
				pos3 = -pos3;
				conj = true;
			}

			const tComplex<SrcType> value = conj? dataFS(x,y).conj() : dataFS(x,y);
			const tComplex<SrcType> wgh = weight(x,y);

			{
				const int x0 = std::floor(pos3.x);
				const int y0 = std::floor(pos3.y);
				const int z0 = std::floor(pos3.z);

				for (int dz = 0; dz < 2; dz++)
				for (int dy = 0; dy < 2; dy++)
				for (int dx = 0; dx < 2; dx++)
				{
					const int xg = x0 + dx;
					const int yg = y0 + dy;
					const int zg = z0 + dz;

					if ( xg < wh3
					  && yg >= -h3/2 && yg < h3/2
					  && zg >= -d3/2 && zg < d3/2)
					{
						const int xi = xg;
						const int yi = yg >= 0? yg : yg + h3;
						const int zi = zg >= 0? zg : zg + d3;

						const double fx = 1.0 - std::abs(pos3.x - xg);
						const double fy = 1.0 - std::abs(pos3.y - yg);
						const double fz = 1.0 - std::abs(pos3.z - zg);

						const double m = scale * fx * fy * fz;

						destFS(  xi,yi,zi) += m * value;
						destCTF( xi,yi,zi) += m * wgh;
						destMult(xi,yi,zi) += m;

						if (xi == 0)
						{
							const int yim = (h3 - yi) % h3;
							const int zim = (d3 - zi) % d3;

							destFS(  0,yim,zim) += m * value.conj();
							destCTF( 0,yim,zim) += m * wgh;
							destMult(0,yim,zim) += m;
						}
					}
				}
			}
		}
	}
}


template <typename SrcType, typename DestType, class PointInsertion>
void FourierBackprojection::backprojectSphere_forward(
				const PointInsertion& pointInsertion,
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				double radius,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;
	
	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSphere_forward: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	const gravis::d3Vector u(proj(0,0), proj(0,1), proj(0,2));
	const gravis::d3Vector v(proj(1,0), proj(1,1), proj(1,2));
	const gravis::d3Vector h(proj(2,0), proj(2,1), proj(2,2));
	
	const double D = 1.0 / (2.0 * radius);
	
	for (long int y = 0; y < h2;  y++)
	for (long int x = (y > 0 && y < h2/2? 0 : 1); x < wh2; x++)
	{
		const double xx = x;
		const double yy = y < h2/2? y : y - h2;
		
		gravis::d3Vector pos3 = xx * u + yy * v + D * (xx*xx + yy*yy) * h;		
		
		bool conj = false;
		
		if (pos3.x < 0)
		{
			pos3 = -pos3;			
			conj = true;
		}		
		
		const tComplex<SrcType> value = conj? dataFS(x,y).conj() : dataFS(x,y);
		const SrcType wgh = weight(x,y);
		
		pointInsertion.insert(value, wgh, pos3, destFS, destCTF);
	}
}

template <typename SrcType, typename DestType, class PointInsertion>
void FourierBackprojection::backprojectSphere_dualContrast_forward(
        const PointInsertion& pointInsertion,
        const RawImage<tComplex<SrcType>>& sin_gamma_data,
		const RawImage<tComplex<SrcType>>& cos_gamma_data,
		const RawImage<SrcType>& sin2_weight,
		const RawImage<SrcType>& sin_cos_weight,
		const RawImage<SrcType>& cos2_weight,
		const gravis::d4Matrix& proj,
		double radius,
        bool conjugate,
		RawImage<DualContrastVoxel<DestType>>& dest)
{
	const int wh2 = sin_gamma_data.xdim;
	const int h2  = sin_gamma_data.ydim;
	
	const gravis::d3Vector u(proj(0,0), proj(0,1), proj(0,2));
	const gravis::d3Vector v(proj(1,0), proj(1,1), proj(1,2));
	const gravis::d3Vector h(proj(2,0), proj(2,1), proj(2,2));
	
	const double D = 1.0 / (2.0 * radius);
	
	for (long int y = 0; y < h2;  y++)
	for (long int x = (y > 0 && y < h2/2? 0 : 1); x < wh2; x++)
	{
		const double xx = x;
		const double yy = y < h2/2? y : y - h2;
		
		gravis::d3Vector pos3 = xx * u + yy * v + D * (xx*xx + yy*yy) * h;		
		
		bool conj = conjugate;
		
		if (pos3.x < 0)
		{
			pos3 = -pos3;			
			conj = !conj;
		}		
		
		gravis::t2Vector<tComplex<SrcType>> value;
		value[0] = conj? sin_gamma_data(x,y).conj() : sin_gamma_data(x,y);
		value[1] = conj? cos_gamma_data(x,y).conj() : cos_gamma_data(x,y);
		
		gravis::t3Vector<SrcType> weight;
		weight[0] = sin2_weight(x,y);
		weight[1] = sin_cos_weight(x,y);
		weight[2] = cos2_weight(x,y);
		
		pointInsertion.insert_dualContrast(value, weight, pos3, dest);
	}
}



template <typename SrcType, typename DestType>
void FourierBackprojection::backprojectSlice_forward_wrap(
				const RawImage<tComplex<SrcType>>& dataFS,
				const RawImage<SrcType>& weight,
				const gravis::d4Matrix& proj,
				RawImage<tComplex<DestType>>& destFS,
				RawImage<DestType>& destCTF)
{
	const int wh2 = dataFS.xdim;
	const int h2 = dataFS.ydim;
	
	const int wh3 = destFS.xdim;
	const int h3 = destFS.ydim;
	const int d3 = destFS.zdim;
	
	if (!destCTF.hasSize(wh3, h3, d3))
	{
		REPORT_ERROR_STR("FourierBackprojection::backprojectSlice_forward_wrap: destCTF has wrong size ("
						 << destCTF.getSizeString() << " instead of " << destCTF.getSizeString() << ")");
	}
	
	const gravis::d3Vector u(proj(0,0), proj(0,1), proj(0,2));
	const gravis::d3Vector v(proj(1,0), proj(1,1), proj(1,2));
	
	for (long int y = 0; y < h2;  y++)
	for (long int x = ((y > 0 && y < h2/2)? 0 : 1); x < wh2; x++)
	{
		const double xx = x;
		const double yy = y < h2/2? y : y - h3;
		
		gravis::d3Vector pos3 = xx * u + yy * v;
		
		bool conj = false;
		
		if (pos3.x < 0)
		{
			pos3 = -pos3;
			conj = true;
		}		
		
		const tComplex<SrcType> value = conj? dataFS(x,y).conj() : dataFS(x,y);
		const tComplex<SrcType> wgh = weight(x,y);
		
		const int x0 = std::floor(pos3.x);
		const int y0 = std::floor(pos3.y);
		const int z0 = std::floor(pos3.z);
		
		for (int dz = 0; dz < 2; dz++)
		for (int dy = 0; dy < 2; dy++)
		for (int dx = 0; dx < 2; dx++)
		{	
			const int xg = x0 + dx;			
			const int yg = y0 + dy;
			const int zg = z0 + dz;
			
			int xi = xg;			
			int yi = yg >= 0? yg : yg + h3;
			int zi = zg >= 0? zg : zg + d3;
			
			tComplex<SrcType> value_wrapped = value;
			
			if (xi > wh3)
			{
				xi = 2*wh3 - xi;
				yi = (h3 - yi) % h3;
				zi = (d3 - zi) % d3;
				
				value_wrapped.imag *= -1;
			}
		
			if ( xi >= 0 && xi < wh3
			  && yi >= 0 && yi < h3
			  && zi >= 0 && zi < d3)
			{		
				const double fx = 1.0 - std::abs(pos3.x - xg);
				const double fy = 1.0 - std::abs(pos3.y - yg);
				const double fz = 1.0 - std::abs(pos3.z - zg);
				
				const double m = fx * fy * fz;
				
				destFS( xi,yi,zi) += m * value_wrapped;
				destCTF(xi,yi,zi) += m * wgh;
			}
		}
	}
}

#endif
