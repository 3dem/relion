#ifndef SKELETONIZATION_H
#define SKELETONIZATION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/image/structure_tensor.h>
#include <src/jaz/util/zio.h>
#include <omp.h>

class Skeletonization
{
	public:
		
		typedef enum
		{
			Point = 0,
			Curve = 1,
			Surface = 2
		}
		Mode;
		
		template <typename T>
		static BufferedImage<T> apply(
			const RawImage<T>& img,
			const RawImage<Tensor3x3<float>>& kernel,				
			int num_threads,
			int r = 1);
		
		static BufferedImage<Tensor3x3<float>> getKernel(
			const RawImage<Tensor3x3<float>>& J,
			Mode mode, float tol = 1.f);
		
		static std::vector<gravis::d3Vector> discretize(
			const RawImage<float>& skeleton,
			float minVal);
};

template <typename T>
BufferedImage<T> Skeletonization :: apply(
	const RawImage<T>& img,
	const RawImage<Tensor3x3<float>>& kernel,
	int num_threads,
	int r)
{
	BufferedImage<T> out(img.xdim, img.ydim, img.zdim);
	
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		gravis::f3Matrix K0(kernel(x,y,z).toMatrix());
		
		const T v0 = img(x,y,z);
		
		bool allSmaller = true;
				
		for (int dz = -r; dz <= r; dz++)
		for (int dy = -r; dy <= r; dy++)
		for (int dx = -r; dx <= r; dx++)
		{
			int xx = x + dx;
			int yy = y + dy;
			int zz = z + dz;
			
			if (xx < 0 || xx >= w || yy < 0 || yy >= h || zz < 0 || zz >= d)
			{
				continue;
			}
			
			gravis::f3Vector p(dx,dy,dz);
			
			const double g0 = p.dot(K0*p);
			
			gravis::f3Matrix K1(kernel(xx,yy,zz).toMatrix());
			
			const double g1 = p.dot(K1*p);
			
			if (g0 > 0 && g1 > 0 && img(xx,yy,zz) > v0)
			{
				allSmaller = false;
				break;
			}
		}
		
		out(x,y,z) = allSmaller? v0 : 0;
	}
	
	return out;
}


#endif
