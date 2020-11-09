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

#ifndef STRUCTURE_TENSOR_H
#define STRUCTURE_TENSOR_H

#include "raw_image.h"
#include "gradient.h"
#include "filter.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/math/tensor3x3.h>
#include <src/jaz/math/tensor2x2.h>
#include <src/jaz/util/new_vtk_helper.h>


class StructureTensor
{
    public:
		
		template <class T>
        static BufferedImage<Tensor3x3<T>> compute3D(
				const RawImage<T>& src, double rho, double sigma0, double taper);
		
		template <class T>
        static BufferedImage<Tensor2x2<T>> compute2D(
				const RawImage<T>& src, double rho, double sigma0, double taper);
		
		template <class T>
        static BufferedImage<Tensor3x3<T>> computeEdgeTensor3D(const RawImage<T>& src, double taper);
		
		template <class T>
        static BufferedImage<Tensor2x2<T>> computeEdgeTensor2D(const RawImage<T>& src, double taper);
		
		template <class T>
        static BufferedImage<gravis::t3Vector<T>> computeEigenvalues3D(const RawImage<Tensor3x3<T>>& J);
		
		template <class T>
        static BufferedImage<gravis::t2Vector<T>> computeEigenvalues2D(const RawImage<Tensor2x2<T>>& J);
		
		template <class T>
		static std::vector<BufferedImage<T>> split(const RawImage<Tensor3x3<T>>& J);
		
		template <class T>
		static std::vector<BufferedImage<T>> split(const RawImage<Tensor2x2<T>>& J);
		
		template <class T>
		static std::vector<BufferedImage<T>> split(const RawImage<gravis::t3Vector<T>>& evals);
		
		template <class T>
		static std::vector<BufferedImage<T>> split(const RawImage<gravis::t2Vector<T>>& evals);
		
		template <class T>
		static BufferedImage<Tensor3x3<T>> join3D(const std::vector<BufferedImage<T>>& Js);
		
		template <class T>
		static BufferedImage<Tensor2x2<T>> join2D(const std::vector<BufferedImage<T>>& Js);
		
		template <class T>
		static BufferedImage<Tensor3x3<T>> downsample(BufferedImage<Tensor3x3<T>>& J, double factor);
		
		template <class T>
		static BufferedImage<Tensor2x2<T>> downsample(BufferedImage<Tensor2x2<T>>& J, double factor);
		
		
		static BufferedImage<Tensor3x3<float>> computeNonLinear(
				const RawImage<float>& src, 
				double rho, 
				double sigma0, 
				double taper, 
				int diffusion_iterations, 
				int outer_iterations, 
				bool planar,
				double lambda,
				double stepSize,
				int num_threads);
		
		template <class T>
        static BufferedImage<Tensor3x3<T>> forwardAverage(
				const RawImage<Tensor3x3<T>>& J);
		
		template <class T>
        static BufferedImage<Tensor2x2<T>> forwardAverage(
				const RawImage<Tensor2x2<T>>& J);
		

};

template <class T>
BufferedImage<Tensor3x3<T>> StructureTensor :: compute3D(const RawImage<T>& src, double rho, double sigma0, double taper)
{
	BufferedImage<T> smoothed = ImageFilter::Gauss3D(src, sigma0, true);
    BufferedImage<Tensor3x3<T>> E = computeEdgeTensor3D(smoothed, taper);
	
	std::vector<BufferedImage<T>> Es = split(E);
	
	if (rho > 0)
	{
		for (int i = 0; i < 6; i++)
		{
			Es[i] = ImageFilter::separableGauss3D(Es[i], rho);
		}
	}
	
	return join3D(Es);
}

template <class T>
BufferedImage<Tensor2x2<T>> StructureTensor :: compute2D(const RawImage<T>& src, double rho, double sigma0, double taper)
{
	BufferedImage<T> smoothed = ImageFilter::GaussStack(src, sigma0, true);
    BufferedImage<Tensor2x2<T>> E = computeEdgeTensor2D(smoothed, taper);
	
	std::vector<BufferedImage<T>> Es = split(E);
	
	if (rho > 0)
	{
		for (int i = 0; i < 3; i++)
		{
			Es[i] = ImageFilter::GaussStack(Es[i], rho,true);
		}
	}
	
	return join2D(Es);
}

template <class T>
BufferedImage<Tensor3x3<T>> StructureTensor :: computeEdgeTensor3D(const RawImage<T>& src, double taper)
{
	const int w = src.xdim;
	const int h = src.ydim;
	const int d = src.zdim;
	
    BufferedImage<Tensor3x3<T>> out(w,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
    {
		gravis::t3Vector<T> g = Gradient::sobelGrid3D(src, x, y, z);
		
		T wgh;		
		if (taper > 0.0) wgh = Tapering::getTaperWeight3D(x,y,z,w,h,d,taper);
		else wgh = 1.0;
			
        out(x,y,z) = Tensor3x3<T>::autoDyadicProduct(wgh * g);
    }
	
	return out;
}

template <class T>
BufferedImage<Tensor2x2<T>> StructureTensor :: computeEdgeTensor2D(const RawImage<T>& src, double taper)
{
	const int w = src.xdim;
	const int h = src.ydim;
	const int fc = src.zdim;
	
    BufferedImage<Tensor2x2<T>> out(w,h,fc);
	
	for (int f = 0; f < fc; f++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
    {
		gravis::t2Vector<T> g = Gradient::sobelGrid2D(src, x, y, f);
		
		T wgh;		
		if (taper > 0.0) wgh = Tapering::getTaperWeight2D(x,y,w,h,taper,0.0);
		else wgh = 1.0;
			
        out(x,y,f) = Tensor2x2<T>::autoDyadicProduct(wgh * g);
    }
	
	return out;
}

template <class T>
BufferedImage<gravis::t3Vector<T>> StructureTensor :: computeEigenvalues3D(const RawImage<Tensor3x3<T>>& J)
{
	BufferedImage<gravis::t3Vector<T>> out(J.xdim, J.ydim, J.zdim);
	
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
    gravis::t3Matrix<T> Q;
    gravis::t3Vector<T> v;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
    {
        J(x,y,z).diagonalize(v, Q);
        out(x,y,z) = v;
    }
	
	return out;
}

template <class T>
BufferedImage<gravis::t2Vector<T>> StructureTensor :: computeEigenvalues2D(const RawImage<Tensor2x2<T>>& J)
{
	BufferedImage<gravis::t2Vector<T>> out(J.xdim, J.ydim, J.zdim);
	
	const int w  = J.xdim;
	const int h  = J.ydim;
	const int fc = J.zdim;
	
    gravis::t2Matrix<T> Q;
    gravis::t2Vector<T> v;
	
	for (int f = 0; f < fc; f++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < w;  x++)
    {
        J(x,y,f).diagonalize(v, Q);
        out(x,y,f) = v;
    }
	
	return out;
}

template <class T>
std::vector<BufferedImage<T>> StructureTensor :: split(const RawImage<Tensor3x3<T>>& J)
{
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	std::vector<BufferedImage<T>> out(6);
	
	for (int i = 0; i < 6; i++)
	{
		out[i] = BufferedImage<T>(w,h,d);
	}
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out[0](x,y,z) = J(x,y,z).xx;
		out[1](x,y,z) = J(x,y,z).xy;
		out[2](x,y,z) = J(x,y,z).xz;
		out[3](x,y,z) = J(x,y,z).yy;
		out[4](x,y,z) = J(x,y,z).yz;
		out[5](x,y,z) = J(x,y,z).zz;
	}
	
	return out;
}

template <class T>
std::vector<BufferedImage<T>> StructureTensor :: split(const RawImage<Tensor2x2<T>>& J)
{
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	std::vector<BufferedImage<T>> out(3);
	
	for (int i = 0; i < 3; i++)
	{
		out[i] = BufferedImage<T>(w,h,d);
	}
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out[0](x,y,z) = J(x,y,z).xx;
		out[1](x,y,z) = J(x,y,z).xy;
		out[2](x,y,z) = J(x,y,z).yy;
	}
	
	return out;
}


template <class T>
std::vector<BufferedImage<T>> StructureTensor :: split(const RawImage<gravis::t3Vector<T>>& evals)
{
	const int w = evals.xdim;
	const int h = evals.ydim;
	const int d = evals.zdim;
	
	std::vector<BufferedImage<T>> out(3);
	
	for (int i = 0; i < 3; i++)
	{
		out[i] = BufferedImage<T>(w,h,d);
	}
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out[0](x,y,z) = evals(x,y,z).x;
		out[1](x,y,z) = evals(x,y,z).y;
		out[2](x,y,z) = evals(x,y,z).z;
	}
	
	return out;
}

template <class T>
std::vector<BufferedImage<T>> StructureTensor :: split(const RawImage<gravis::t2Vector<T>>& evals)
{
	const int w = evals.xdim;
	const int h = evals.ydim;
	const int d = evals.zdim;
	
	std::vector<BufferedImage<T>> out(2);
	
	for (int i = 0; i < 2; i++)
	{
		out[i] = BufferedImage<T>(w,h,d);
	}
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out[0](x,y,z) = evals(x,y,z).x;
		out[1](x,y,z) = evals(x,y,z).y;
	}
	
	return out;
}

template <class T>
BufferedImage<Tensor3x3<T>> StructureTensor :: join3D(const std::vector<BufferedImage<T>>& Js)
{
	const int w = Js[0].xdim;
	const int h = Js[0].ydim;
	const int d = Js[0].zdim;
	
	BufferedImage<Tensor3x3<T>> out(w,h,d);
				
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y,z).xx = Js[0](x,y,z);
		out(x,y,z).xy = Js[1](x,y,z);
		out(x,y,z).xz = Js[2](x,y,z);
		out(x,y,z).yy = Js[3](x,y,z);
		out(x,y,z).yz = Js[4](x,y,z);
		out(x,y,z).zz = Js[5](x,y,z);
	}
	
	return out;
}

template <class T>
BufferedImage<Tensor2x2<T>> StructureTensor :: join2D(const std::vector<BufferedImage<T>>& Js)
{
	const int w = Js[0].xdim;
	const int h = Js[0].ydim;
	const int d = Js[0].zdim;
	
	BufferedImage<Tensor2x2<T>> out(w,h,d);
				
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y,z).xx = Js[0](x,y,z);
		out(x,y,z).xy = Js[1](x,y,z);
		out(x,y,z).yy = Js[2](x,y,z);
	}
	
	return out;
}

template <class T>
BufferedImage<Tensor3x3<T>> StructureTensor :: downsample(BufferedImage<Tensor3x3<T>>& J, double factor)
{
	const int w0 = J.xdim;
	const int h0 = J.ydim;
	const int d0 = J.zdim;
	
	const int w1 = (int) (w0 / factor);
	const int h1 = (int) (h0 / factor);
	const int d1 = (int) (d0 / factor);
	
	std::vector<BufferedImage<T>> components = StructureTensor::split(J);
	
	for (int i = 0; i < components.size(); i++)
	{
		//components[i] = Resampling::downsampleGaussFilt_3D_full(components[i], w1, h1, d1, factor);
		components[i] = Resampling::subsample_3D_full(components[i], w1, h1, d1, factor);
	}
	
	return join3D(components);
}

template <class T>
BufferedImage<Tensor2x2<T>> StructureTensor :: downsample(BufferedImage<Tensor2x2<T>>& J, double factor)
{
	const int w0 = J.xdim;
	const int h0 = J.ydim;
	const int d0 = J.zdim;
	
	const int w1 = (int) (w0 / factor);
	const int h1 = (int) (h0 / factor);
	
	std::vector<BufferedImage<T>> components = StructureTensor::split(J);
	
	for (int i = 0; i < components.size(); i++)
	{
		//components[i] = Resampling::downsampleGaussFilt_3D_full(components[i], w1, h1, d1, factor);
		components[i] = Resampling::subsample_2D_full(components[i], w1, h1, factor);
	}
	
	return join2D(components);
}


template <class T>
BufferedImage<Tensor3x3<T>> StructureTensor :: forwardAverage(
		const RawImage<Tensor3x3<T>>& J)
{
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	BufferedImage<Tensor3x3<T>> out(w,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		Tensor3x3<T> sum(0);
		T wgh(0);
		
		for (int xx = x; xx <= x+1 && xx < w; xx++)
		for (int yy = y; yy <= y+1 && yy < h; yy++)
		for (int zz = z; zz <= z+1 && zz < d; zz++)
		{
			sum += J(xx,yy,zz);
			wgh++;
		}
		
		out(x,y,z) = sum / wgh;
	}
	
	return out;
}

template <class T>
BufferedImage<Tensor2x2<T>> StructureTensor :: forwardAverage(
		const RawImage<Tensor2x2<T>>& J)
{
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	BufferedImage<Tensor2x2<T>> out(w,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		Tensor2x2<T> sum(0);
		T wgh(0);
		
		for (int xx = x; xx <= x+1 && xx < w; xx++)
		for (int yy = y; yy <= y+1 && yy < h; yy++)
		{
			sum += J(xx,yy,z);
			wgh++;
		}
		
		out(x,y,z) = sum / wgh;
	}
	
	return out;
}

#endif
