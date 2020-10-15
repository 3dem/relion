#ifndef IMAGE_GRADIENT_H
#define IMAGE_GRADIENT_H

#include "buffered_image.h"
#include <src/jaz/gravis/t3Vector.h>

class Gradient
{
	public:
		
		template<typename T> 
		static inline gravis::t3Vector<T> sobelGrid3D(const RawImage<T>& img, int x, int y, int z);
		
		template<typename T> 
		static inline gravis::t2Vector<T> sobelGrid2D(const RawImage<T>& img, int x, int y, int z);
		
		template<typename T> 
		static inline gravis::t3Vector<T> central3D(const RawImage<T>& img, int x, int y, int z);
		
		template<typename T>
		static inline void forward3D_inSitu(const RawImage<T>& img, RawImage<gravis::t3Vector<T>>& dest);
		
		
		template<typename T>
		static BufferedImage<gravis::t3Vector<T>> computeGradImage(const RawImage<T>& img);
		
		template<typename T>
		static BufferedImage<T> computeGradL2Image(const RawImage<T>& img);
		
};

template<typename T> 
inline gravis::t3Vector<T> Gradient::sobelGrid3D(const RawImage<T>& img, int x, int y, int z)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	const int xn = x > 0? x - 1 : 0;
	const int xp = x >= w-1? w - 1 : x + 1;
	
	const int yn = y > 0? y - 1 : 0;
	const int yp = y >= h-1? h - 1 : y + 1;
	
	const int zn = z > 0? z - 1 : 0;
	const int zp = z >= d-1? d - 1 : z + 1;
	
	const T* v0 = img.data;
	
	gravis::t3Vector<size_t> strides(
		1, 
		img.xdim, 
		img.xdim * img.ydim);
	
	gravis::t3Matrix<size_t> indices(
		xn, x, xp,
		yn, y, yp,
		zn, z, zp);
	
	gravis::t3Vector<T> out;
	
	for (int dim = 0; dim < 3; dim++)
	{
		out[dim] = T(0);
		
		const int p = (dim + 1) % 3;
		const int q = (dim + 2) % 3;
		
		for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
		{
			const double wgh = (2 - i*i) * (2 - j*j) / 32.0;
			
			const size_t ind_dl = indices(dim, 0);
			const size_t ind_dr = indices(dim, 2);
			const size_t ind_p = indices(p, i+1);
			const size_t ind_q = indices(q, j+1);
			
			const T vl = v0[ind_p * strides[p] + ind_q * strides[q] + ind_dl * strides[dim]];
			const T vr = v0[ind_p * strides[p] + ind_q * strides[q] + ind_dr * strides[dim]];
					
			out[dim] += wgh * (vr - vl);
		}
	}
		 
	return out;
}


template<typename T> 
inline gravis::t2Vector<T> Gradient::sobelGrid2D(const RawImage<T>& img, int x, int y, int z)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
	const int xn = x > 0? x-1 : 0;
	const int xp = x >= w-1? w-1 : x+1;
	
	const int yn = y > 0? y-1 : 0;
	const int yp = y >= h-1? h-1 : y+1;
		
	return gravis::t2Vector<T>(
		(   (img(xp,yn,z) + T(2) * img(xp,y,z) + img(xp,yp,z))
		  - (img(xn,yn,z) + T(2) * img(xn,y,z) + img(xn,yp,z)) ) / T(8),
		(   (img(xn,yp,z) + T(2) * img(x,yp,z) + img(xp,yp,z))
		  - (img(xn,yn,z) + T(2) * img(x,yn,z) + img(xp,yn,z)) ) / T(8));
}


template<typename T> 
inline gravis::t3Vector<T> Gradient::central3D(const RawImage<T>& img, int x, int y, int z)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	const int xn = x > 0? x - 1 : 0;
	const int xp = x >= w-1? w - 1 : x + 1;
	
	const int yn = y > 0? y - 1 : 0;
	const int yp = y >= h-1? h - 1 : y + 1;
	
	const int zn = z > 0? z - 1 : 0;
	const int zp = z >= d-1? d - 1 : z + 1;
	
	return gravis::t3Vector<T>(
		0.5 * (img(xp,y,z) - img(xn,y,z)),
		0.5 * (img(x,yp,z) - img(x,yn,z)),
		0.5 * (img(x,y,zp) - img(x,y,zn)));
}

template<typename T> 
inline void Gradient::forward3D_inSitu(
		const RawImage<T>& img, RawImage<gravis::t3Vector<T>>& dest)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const int xp = x >= w - 1?  w - 2 : x + 1;
		const int x0 = x >= w - 1?  w - 1 : x;
		
		const int yp = y >= h - 1?  h - 2 : y + 1;
		const int y0 = y >= h - 1?  h - 1 : y;
		
		const int zp = z >= d - 1?  d - 2 : z + 1;
		const int z0 = z >= d - 1?  d - 1 : z;
		
		dest(x,y,z) = gravis::t3Vector<T>(
			(img(xp,y,z) - img(x0,y,z)),
			(img(x,yp,z) - img(x,y0,z)),
			(img(x,y,zp) - img(x,y,z0)));
	}
}

template<typename T>
BufferedImage<gravis::t3Vector<T>> Gradient::computeGradImage(const RawImage<T>& img)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	BufferedImage<gravis::t3Vector<T>> out(w,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y,z) = sobelGrid3D(img, x, y, z);		
	}
	
	return out;
}

template<typename T>
BufferedImage<T> Gradient::computeGradL2Image(const RawImage<T>& img)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	BufferedImage<gravis::t3Vector<T>> grad = computeGradImage(img);
	
	BufferedImage<T> out(w,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y,z) = grad(x,y,z).norm2();		
	}
	
	return out;
}


#endif
