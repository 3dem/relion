/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * C-CINA Basel
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

#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_voxel.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_solution.h>
#include "extraction.h"
#include "tomo_stack.h"


class Reconstruction
{
    public:

		template <typename T>
		static void griddingCorrect3D(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<T>& psfImgFS,
				BufferedImage<T>& out,
				bool center,
				int num_threads);
		
		template <typename T>
		static void griddingCorrect_dualContrast(
				BufferedImage<DualContrastVoxel<T>>& image,
				BufferedImage<T>& psfImgFS,
				BufferedImage<DualContrastVoxel<T>>& out,
				bool center,
				int num_threads);
		
		template <typename T>
		static void griddingCorrectSinc2_dualContrast(
				BufferedImage<DualContrastVoxel<T>>& image,
				BufferedImage<DualContrastVoxel<T>>& out,
				bool center,
				int num_threads);
		
		template <typename T>
		static void griddingCorrect3D_sinc2(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<T>& out,
				bool center,
				int num_threads);

		template <typename T>
		static void ctfCorrect3D_Wiener(
				BufferedImage<T>& img,
				BufferedImage<T>& ctfImgFS,
				BufferedImage<T>& out,
				double WienerOffset,
				int num_threads = 1);

		template <typename T>
		static void ctfCorrect3D_heuristic(
				BufferedImage<T>& img,
				BufferedImage<T>& ctfImgFS,
				BufferedImage<T>& out,
				double weightFraction = 0.001,
				int num_threads = 1);

		template <typename T>
		static DualContrastSolution<T> solveDualContrast(
				BufferedImage<DualContrastVoxel<T>>& image,
				double SNR = 1,
				double lambda = 0,
				bool isotropicWiener = true,
				int num_threads = 1);
		
		template <typename T>
		static void taper(
				BufferedImage<T>& img,
				double dist,
				bool do_subtract_mean,
				int num_threads);
		
		template <typename T>
		static void GaussEnvelope(
				BufferedImage<T>& img,
				double sigma,
				bool do_subtract_mean,
				int num_threads);
		
		
		template <typename T>
		static void correct3D_RS(
				BufferedImage<T>& dataImgRS,
				BufferedImage<T>& psfImgRS,
				BufferedImage<T>& out,
				double WienerOffset,
				int num_threads);
		
		// used after a forward projection:
		template <typename T>
		static void correctStack(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<tComplex<T>>& psfImgFS,
				RawImage<T>& out,
				bool center, 
				int num_threads);
		
		template <typename T>
		static void correctStackFS(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<tComplex<T>>& psfImgFS,
				BufferedImage<tComplex<T>>& out,
				bool center, 
				int num_threads);
		
		// used after a cross projection:
		template <typename T>
		static void correctStack(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<T>& psfImgFS,
				BufferedImage<T>& ctfImgFS,
				RawImage<T>& out,
				double WienerOffset,
				bool center, 
				int num_threads);

};

template <typename T>
void Reconstruction :: griddingCorrect3D(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<T>& psfImgFS,
		BufferedImage<T>& out,
		bool center, 
		int num_threads)
{
	const long int wh = dataImgFS.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = dataImgFS.ydim;
	const long int d = dataImgFS.zdim;
	
	if (center)
	{
		#pragma omp parallel for num_threads(num_threads)
		for (long int z = 0; z < d;  z++)
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			dataImgFS(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}
	
	// gridding correction
	
	BufferedImage<T> dataImg(w,h,d);	
	FFT::inverseFourierTransform(dataImgFS, dataImg, FFT::Both);

	BufferedImage<tComplex<T>> psfImgFS_complex(wh,h,d);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		psfImgFS_complex(x,y,z) = psfImgFS(x,y,z);
		
		if (center) 
		{
			psfImgFS_complex(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}
	
	BufferedImage<T> dataEnv(w,h,d);	
	FFT::inverseFourierTransform(psfImgFS_complex, dataEnv, FFT::Both);
		
	const double eps = 1e-10;
	
	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		if (dataEnv(x,y,z) > eps)
		{
			out(x,y,z) = dataImg(x,y,z) / dataEnv(x,y,z);
		}
	}
}

template <typename T>
void Reconstruction :: griddingCorrect_dualContrast(
		BufferedImage<DualContrastVoxel<T>>& image,
		BufferedImage<T>& psfImgFS,
		BufferedImage<DualContrastVoxel<T>>& out,
		bool center,
		int num_threads)
{
	const long int wh = image.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = image.ydim;
	const long int d = image.zdim;

	if (!out.hasEqualSize(image))
	{
		REPORT_ERROR_STR("Reconstruction :: griddingCorrect_dualContrast: unequal input and output image sizes");
	}

	if (center)
	{
		#pragma omp parallel for num_threads(num_threads)
		for (long int z = 0; z < d;  z++)
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			image(x,y,z).data_sin *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
			image(x,y,z).data_cos *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}

	BufferedImage<tComplex<T>> sinImgFS(wh,h,d), cosImgFS(wh,h,d);

	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		sinImgFS(x,y,z) = image(x,y,z).data_sin;
		cosImgFS(x,y,z) = image(x,y,z).data_cos;
	}

	BufferedImage<T> sinImgRS(w,h,d), cosImgRS(w,h,d);

	FFT::inverseFourierTransform(sinImgFS, sinImgRS, FFT::Both);
	FFT::inverseFourierTransform(cosImgFS, cosImgRS, FFT::Both);

	BufferedImage<tComplex<T>> psfImgFS_complex(wh,h,d);

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		psfImgFS_complex(x,y,z) = psfImgFS(x,y,z);

		if (center)
		{
			psfImgFS_complex(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}

	BufferedImage<T> dataEnv(w,h,d);
	FFT::inverseFourierTransform(psfImgFS_complex, dataEnv, FFT::Both);

	const double eps = 1e-10;

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		if (dataEnv(x,y,z) > eps)
		{
			sinImgRS(x,y,z) /= dataEnv(x,y,z);
			cosImgRS(x,y,z) /= dataEnv(x,y,z);
		}
	}

	FFT::FourierTransform(sinImgRS, sinImgFS, FFT::Both);
	FFT::FourierTransform(cosImgRS, cosImgFS, FFT::Both);

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		out(x,y,z) = image(x,y,z);
		out(x,y,z).data_sin = sinImgFS(x,y,z);
		out(x,y,z).data_cos = cosImgFS(x,y,z);
	}
}

template <typename T>
void Reconstruction :: griddingCorrectSinc2_dualContrast(
		BufferedImage<DualContrastVoxel<T>>& image,
		BufferedImage<DualContrastVoxel<T>>& out,
		bool center,
		int num_threads)
{
	const long int wh = image.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = image.ydim;
	const long int d = image.zdim;

	if (!out.hasEqualSize(image))
	{
		REPORT_ERROR_STR("Reconstruction :: griddingCorrectSinc2_dualContrast: unequal input and output image sizes");
	}

	if (center)
	{
		#pragma omp parallel for num_threads(num_threads)
		for (long int z = 0; z < d;  z++)
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			image(x,y,z).data_sin *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
			image(x,y,z).data_cos *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}

	BufferedImage<tComplex<T>> sinImgFS(wh,h,d), cosImgFS(wh,h,d);

	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		sinImgFS(x,y,z) = image(x,y,z).data_sin;
		cosImgFS(x,y,z) = image(x,y,z).data_cos;
	}

	BufferedImage<T> sinImgRS(w,h,d), cosImgRS(w,h,d);

	FFT::inverseFourierTransform(sinImgFS, sinImgRS, FFT::Both);
	FFT::inverseFourierTransform(cosImgFS, cosImgRS, FFT::Both);


	const double eps = 1e-10;

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		const double xx = x - w/2;
		const double yy = y - h/2;
		const double zz = z - d/2;
		
		if (xx == 0 && yy == 0 && zz == 0)
		{
			// sinc at 0 is 1
		}
		else
		{		
			const double r = sqrt(xx*xx + yy*yy + zz*zz);
			const double d = r / w;
			const double sinc = sin(PI * d) / (PI * d);
			const double sinc2 = sinc * sinc;
			
			if (sinc2 > eps)
			{
				sinImgRS(x,y,z) /= sinc2;
				cosImgRS(x,y,z) /= sinc2;
			}
			else
			{
				sinImgRS(x,y,z) /= eps;
				cosImgRS(x,y,z) /= eps;
			}			
		}
	}

	FFT::FourierTransform(sinImgRS, sinImgFS, FFT::Both);
	FFT::FourierTransform(cosImgRS, cosImgFS, FFT::Both);

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		out(x,y,z) = image(x,y,z);
		out(x,y,z).data_sin = sinImgFS(x,y,z);
		out(x,y,z).data_cos = cosImgFS(x,y,z);
	}
}


template <typename T>
void Reconstruction :: griddingCorrect3D_sinc2(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<T>& out,
		bool center, 
		int num_threads)
{
	const long int wh = dataImgFS.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = dataImgFS.ydim;
	const long int d = dataImgFS.zdim;
	
	if (center)
	{
		#pragma omp parallel for num_threads(num_threads)
		for (long int z = 0; z < d;  z++)
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			dataImgFS(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}
	
	// gridding correction
	
	BufferedImage<T> dataImg(w,h,d);	
	FFT::inverseFourierTransform(dataImgFS, dataImg, FFT::Both);
	
	const double eps = 1e-2;
	
	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		const double xx = x - w/2;
		const double yy = y - h/2;
		const double zz = z - d/2;
		
		if (xx == 0 && yy == 0 && zz == 0)
		{
			out(x,y,z) = dataImg(x,y,z);
		}
		else
		{		
			const double r = sqrt(xx*xx + yy*yy + zz*zz);
			const double d = r / w;
			const double sinc = sin(PI * d) / (PI * d);
			const double sinc2 = sinc * sinc;

			if (sinc2 < eps || d > 1.0)
			{
				out(x,y,z) = dataImg(x,y,z) / eps;
			}
			else
			{
				out(x,y,z) = dataImg(x,y,z) / sinc2;
			}
		}
	}
}

template <typename T>
void Reconstruction :: ctfCorrect3D_Wiener(
		BufferedImage<T>& img,
		BufferedImage<T>& ctfImgFS,
		BufferedImage<T>& out,
		double WienerOffset,
		int num_threads)
{
	const long int w = img.xdim;
	const long int wh = w/2 + 1;
	const long int h = img.ydim;
	const long int d = img.zdim;

	BufferedImage<tComplex<T>> dataImgCorrFS(wh,h,d);
	FFT::FourierTransform(img, dataImgCorrFS, FFT::Both);

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		dataImgCorrFS(x,y,z) /= ctfImgFS(x,y,z) + WienerOffset;
	}

	FFT::inverseFourierTransform(dataImgCorrFS, out, FFT::Both);
}

template <typename T>
void Reconstruction :: ctfCorrect3D_heuristic(
		BufferedImage<T>& img,
		BufferedImage<T>& ctfImgFS,
		BufferedImage<T>& out,
		double weightFraction,
		int num_threads)
{
	const long int w = img.xdim;
	const long int wh = w/2 + 1;
	const long int h = img.ydim;
	const long int d = img.zdim;

	BufferedImage<tComplex<T>> dataImgCorrFS(wh,h,d);
	FFT::FourierTransform(img, dataImgCorrFS, FFT::Both);

	std::vector<T> rad_avg = RadialAvg::fftwHalf_3D_lin(ctfImgFS);

	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		/*const T avgWeight = RadialAvg::interpolate_FftwHalf_3D_lin(
					x, y, z, w, h, d, rad_avg);*/

		const int s = rad_avg.size();

		const double xx = x;
		const double yy = (y >= h/2)? y - h: y;
		const double zz = (z >= d/2)? z - d: z;

		const double r = sqrt(xx * xx + yy * yy + zz * zz);

		const int r0 = (int) r;
		const int r1 = r0 + 1;

		if (r1 >= s)
		{
			dataImgCorrFS(x,y,z) = T(0);
		}
		else
		{
			const T v1 = rad_avg[r1];
			const T v0 = rad_avg[r0];

			const double f = r - r0;

			const T avgWeight = f * v1 + (1.0 - f) * v0;

			const T threshold = avgWeight * weightFraction;

			T wg = ctfImgFS(x,y,z);

			if (wg < threshold)
			{
				wg = threshold;
			}

			if (wg > 0.0)
			{
				dataImgCorrFS(x,y,z) /= wg;
			}
		}
	}

	FFT::inverseFourierTransform(dataImgCorrFS, out, FFT::Both);
}

template<typename T>
DualContrastSolution<T>
	Reconstruction::solveDualContrast(
		BufferedImage<DualContrastVoxel<T>>& image,
		double SNR,
		double lambda,
		bool isotropicWiener,
		int num_threads)
{
	const long int wh = image.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = image.ydim;
	const long int d = image.zdim;

	const double WienerOffset = SNR > 0.0? 1.0 / SNR : 0.0;

	BufferedImage<tComplex<T>> phaseOutFS(wh,h,d), ampOutFS(wh,h,d);
	BufferedImage<T> conditionNumber(wh,h,d);


	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		typename DualContrastVoxel<T>::Solution X = image(x,y,z).solve(
					WienerOffset, lambda, isotropicWiener);

		phaseOutFS(x,y,z) = X.phase;
		ampOutFS(x,y,z)   = X.amplitude;

		conditionNumber(x,y,z) = X.conditionNumber;
	}

	DualContrastSolution<T> out(w,h,d);

	FFT::inverseFourierTransform(phaseOutFS, out.phase, FFT::Both);
	FFT::inverseFourierTransform(ampOutFS, out.amplitude, FFT::Both);


	std::vector<int> shellVolume(wh, 0.0);

	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		const double r = RadialAvg::get1DIndex(x,y,z, w,h,d);
		const int ri = (int)(r + 0.5);

		if (ri >= wh) continue;

		const double c = conditionNumber(x,y,z);

		if (out.conditionPerShell[ri].maximum < c)
		{
			out.conditionPerShell[ri].maximum = c;
		}

		if (out.conditionPerShell[ri].minimum > c)
		{
			out.conditionPerShell[ri].minimum = c;
		}

		out.conditionPerShell[ri].mean += c;

		shellVolume[ri]++;
	}

	for (int r = 0; r < wh; r++)
	{
		if (shellVolume[r] > 0)
		{
			out.conditionPerShell[r].mean /= shellVolume[r];
		}
	}

	std::vector<double> shellVariance(wh, 0.0);

	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		const double r = RadialAvg::get1DIndex(x,y,z, w,h,d);
		const int ri = (int)(r + 0.5);

		if (ri >= wh) continue;

		const double d = conditionNumber(x,y,z) - out.conditionPerShell[ri].mean;

		shellVariance[ri] += d * d;
	}

	for (int r = 0; r < wh; r++)
	{
		if (shellVolume[r] > 1)
		{
			shellVariance[r] /= (shellVolume[r] - 1);
			out.conditionPerShell[r].std_deviation = sqrt(shellVariance[r]);
		}
	}

	return out;
}

template <typename T>
void Reconstruction :: taper(
		BufferedImage<T>& img,
		double dist,
		bool do_subtract_mean,
		int num_threads)
{
	double innerAvg = 0.0;
	double innerWgh = 0.0;
	
	const int s = img.xdim;
	const double m = s/2.0;
	
	if (do_subtract_mean)
	{
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double dx = x - m;
			const double dy = y - m;
			const double dz = z - m;
			
			const double r = sqrt(dx*dx + dy*dy + dz*dz);
			
			if (r < m)
			{
				double c = r < m - dist? 1.0 : 0.5 - 0.5 * cos(PI * (m - r) / dist);
				
				innerAvg += c * img(x,y,z);
				innerWgh += c;
			}
		}
		
		innerAvg /= innerWgh;
	}
	else
	{
		innerAvg = 0.0;
	}
	
	
	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double dx = x - m;
		const double dy = y - m;
		const double dz = z - m;
		
		const double r = sqrt(dx*dx + dy*dy + dz*dz);
		
		if (r < m)
		{
			double c = r < m - dist? 1.0 : 0.5 - 0.5 * cos(PI * (m - r) / dist);
			
			img(x,y,z) = c * (img(x,y,z) - innerAvg);
		}
		else
		{
			img(x,y,z) = 0;
		}
	}
}


template <typename T>
void Reconstruction :: GaussEnvelope(
		BufferedImage<T>& img,
		double sigma,
		bool do_subtract_mean,
		int num_threads)
{
	if (sigma < 0.0) return;
	
	double innerAvg = 0.0;
	double innerWgh = 0.0;
	
	const int s = img.xdim;
	const double m = s/2.0;
	const double beta = -1.0 / (2.0 * sigma * sigma);
	
	if (do_subtract_mean)
	{
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double dx = x - m;
			const double dy = y - m;
			const double dz = z - m;
			
			const double c = exp(beta * (dx*dx + dy*dy + dz*dz));
			
			innerAvg += c * img(x,y,z);
			innerWgh += c;
		}
		
		innerAvg /= innerWgh;
	}
	else
	{
		innerAvg = 0.0;
	}
	
	
	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double dx = x - m;
		const double dy = y - m;
		const double dz = z - m;
		
		const double c = exp(beta * (dx*dx + dy*dy + dz*dz));
		
		img(x,y,z) = c * (img(x,y,z) - innerAvg);
	}
}

template <typename T>
void Reconstruction :: correct3D_RS(
		BufferedImage<T>& dataImgRS,
		BufferedImage<T>& psfImgRS,
		BufferedImage<T>& out,
		double WienerOffset,
		int num_threads)
{
	BufferedImage<fComplex> dataFS, psfFS;
	
	FFT::FourierTransform(dataImgRS, dataFS, FFT::Both);
	FFT::FourierTransform(psfImgRS, psfFS, FFT::Both);
		
	#pragma omp parallel for num_threads(num_threads)
	for (int z = 0; z < dataFS.zdim; z++)
	for (int y = 0; y < dataFS.ydim; y++)
	for (int x = 0; x < dataFS.xdim; x++)
	{
		const double wa = psfFS(x,y,z).abs() + WienerOffset;
		
		if (wa > 0.0) dataFS(x,y,z) /= wa;
	}
	
	FFT::inverseFourierTransform(dataFS, out, FFT::Both);
}

template <typename T>
void Reconstruction :: correctStack(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<tComplex<T>>& psfImgFS,
		RawImage<T>& out,
		bool center, 
		int num_threads)
{
	const int wh = dataImgFS.xdim;
	const int w  = 2 * (wh - 1);
	const int h  = dataImgFS.ydim;
	const int fc = dataImgFS.zdim;
	
	BufferedImage<T> data(w,h,fc), psf(w,h,fc);
	
	NewStackHelper::inverseFourierTransformStack(dataImgFS, data, center);
	NewStackHelper::inverseFourierTransformStack(psfImgFS, psf, center);
		
	#pragma omp parallel for num_threads(num_threads)
	for (long int f = 0; f < fc; f++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		out(x,y,f) = data(x,y,f) / psf(x,y,f);
	}
}

template <typename T>
void Reconstruction :: correctStackFS(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<tComplex<T>>& psfImgFS,
		BufferedImage<tComplex<T>>& out,
		bool center, 
		int num_threads)
{
	const int wh = dataImgFS.xdim;
	const int w  = 2 * (wh - 1);
	const int h  = dataImgFS.ydim;
	const int fc = dataImgFS.zdim;
	
	BufferedImage<T> temp(w,h,fc);
	
	correctStack(dataImgFS, psfImgFS, temp, center, num_threads);
	
	out = NewStackHelper::FourierTransformStack(temp, center, num_threads);
}


template <typename T>
void Reconstruction :: correctStack(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<T>& psfImgFS,
		BufferedImage<T>& ctfImgFS,
		RawImage<T>& out,
		double WienerOffset,
		bool center, 
		int num_threads)
{
	const long int wh = dataImgFS.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = dataImgFS.ydim;
	const long int fc = dataImgFS.zdim;
	
	// deconvolve by Fourier-space 2D-PSF:
	
	std::vector<BufferedImage<T>> 
			tempDataRS(num_threads, BufferedImage<T>(w,h)),
			tempEnvRS(num_threads, BufferedImage<T>(w,h));
	
	std::vector<BufferedImage<tComplex<T>>> 
			tempDataFS(num_threads, BufferedImage<tComplex<T>>(wh,h)),
			tempEnvFS(num_threads, BufferedImage<tComplex<T>>(wh,h));
	
	#pragma omp parallel for num_threads(num_threads)			
	for (long int f = 0; f < fc;  f++)
	{
		const int t = omp_get_thread_num();
		
		if (center)
		{
			for (long int y = 0; y < h;  y++)
			for (long int x = 0; x < wh; x++)
			{
				tempDataFS[t](x,y) = dataImgFS(x,y,f) * (1 - 2*(x%2)) * (1 - 2*(y%2));
				tempEnvFS[t](x,y) = psfImgFS(x,y,f) * (1 - 2*(x%2)) * (1 - 2*(y%2));
			}
		}
		else
		{
			for (long int y = 0; y < h;  y++)
			for (long int x = 0; x < wh; x++)
			{
				tempDataFS[t](x,y) = dataImgFS(x,y,f);
				tempEnvFS[t](x,y) = psfImgFS(x,y,f);
			}
		}
		
		FFT::inverseFourierTransform(tempDataFS[t], tempDataRS[t], FFT::Both);
		FFT::inverseFourierTransform(tempEnvFS[t], tempEnvRS[t], FFT::Both);
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			out(x,y,f) = tempDataRS[t](x,y) / tempEnvRS[t](x,y);
		}
	}
	
	// divide by ctf in Fourier space
	
	#pragma omp parallel for num_threads(num_threads)
	for (long int f = 0; f < fc;  f++)
	{
		const int t = omp_get_thread_num();
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			tempDataRS[t](x,y) = out(x,y,f);
		}
		
		FFT::FourierTransform(tempDataRS[t], tempDataFS[t], FFT::Both);
		
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			tempDataFS[t](x,y) /= (ctfImgFS(x,y,f) + WienerOffset);
		}
		
		FFT::inverseFourierTransform(tempDataFS[t], tempDataRS[t], FFT::Both);
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			out(x,y,f) = tempDataRS[t](x,y);
		}
	}

}

#endif
