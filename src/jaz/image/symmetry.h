#ifndef JAZ_SYMMETRY_H
#define JAZ_SYMMETRY_H

#include "buffered_image.h"
#include "interpolation.h"
#include <src/symmetries.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_voxel.h>
#include <omp.h>

class Symmetry
{
	public:
		
		template <typename T>
		static BufferedImage<T> symmetrise_FS_real(
				const RawImage<T>& img, std::string symmGroup, int num_threads);

		template <typename T>
		static BufferedImage<tComplex<T>> symmetrise_FS_complex(
				const RawImage<tComplex<T>>& img, std::string symmGroup, int num_threads);

		template <typename T>
		static BufferedImage<DualContrastVoxel<T>> symmetrise_dualContrast(
				const RawImage<DualContrastVoxel<T>>& img, std::string symmGroup, int num_threads);

		static std::vector<gravis::d4Matrix> getMatrices(
				std::string symmGroup);
};

template <typename T>
BufferedImage<T> Symmetry::symmetrise_FS_real(
		const RawImage<T>& img, std::string symmGroup, int num_threads)
{
	std::vector<gravis::d4Matrix> R = getMatrices(symmGroup);
	
	const int wh = img.xdim;
	const int h  = img.ydim;
	const int d  = img.zdim;
	const int sc = R.size();
	
	BufferedImage<T> out(wh,h,d);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = y < h/2? y : y - h;
		const double zz = z < d/2? z : z - d;
		
		T accum = img(x,y,z);
		
		for (int sym = 0; sym < sc; sym++)
		{
			gravis::d4Vector p = R[sym] * gravis::d4Vector(xx,yy,zz, 0.0);
			
			accum += Interpolation::linearXYZ_FftwHalf_real(img, p.x, p.y, p.z);
		}
		
		out(x,y,z) = accum;														   
	}
	
	return out;
}

template <typename T>
BufferedImage<tComplex<T>> Symmetry::symmetrise_FS_complex(
		const RawImage<tComplex<T>>& img, std::string symmGroup, int num_threads)
{
	std::vector<gravis::d4Matrix> R = getMatrices(symmGroup);

	const int wh = img.xdim;
	const int h  = img.ydim;
	const int d  = img.zdim;
	const int sc = R.size();

	BufferedImage<tComplex<T>> out(wh,h,d);

	#pragma omp parallel for num_threads(num_threads)
	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = y < h/2? y : y - h;
		const double zz = z < d/2? z : z - d;

		tComplex<T> accum = img(x,y,z);

		for (int sym = 0; sym < sc; sym++)
		{
			gravis::d4Vector p = R[sym] * gravis::d4Vector(xx,yy,zz, 0.0);

			accum += Interpolation::linearXYZ_FftwHalf_complex(img, p.x, p.y, p.z);
		}

		out(x,y,z) = accum;
	}

	return out;
}

template <typename T>
BufferedImage<DualContrastVoxel<T>> Symmetry::symmetrise_dualContrast(
		const RawImage<DualContrastVoxel<T>>& img, std::string symmGroup, int num_threads)
{
	std::vector<gravis::d4Matrix> R = getMatrices(symmGroup);

	const int wh = img.xdim;
	const int h  = img.ydim;
	const int d  = img.zdim;
	const int sc = R.size();

	BufferedImage<DualContrastVoxel<T>> out(wh,h,d);

	#pragma omp parallel for num_threads(num_threads)
	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = y < h/2? y : y - h;
		const double zz = z < d/2? z : z - d;

		DualContrastVoxel<T> accum = img(x,y,z);

		for (int sym = 0; sym < sc; sym++)
		{
			gravis::d4Vector p = R[sym] * gravis::d4Vector(xx,yy,zz, 0.0);

			accum += Interpolation::linearXYZ_FftwHalf_generic(img, p.x, p.y, p.z);
		}

		out(x,y,z) = accum;
	}

	return out;
}

#endif

