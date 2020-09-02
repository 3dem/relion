#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/image/structure_tensor.h>
#include <src/jaz/util/zio.h>
#include <omp.h>

class Diffusion
{
	public:
		
		/* Diffuses image data along the diffusion tensor D. */
		template <typename T>
		static BufferedImage<T> diffuse(
			const RawImage<T>& data,
			const RawImage<Tensor3x3<float>>& D,
			float step, int iterations, int num_threads);
		
		/* Diffuses by only moving the function value upward.
		   Useful to compute maxima along a certain direction or pair of directions. */
		template <typename T>
		static BufferedImage<T> diffuseUpward(
			const RawImage<T>& data,
			const RawImage<Tensor3x3<float>>& D,
			float step, int iterations, int num_threads);
		
		/* Diffuses the image while maintaining the value from 'data' where 'weight' is large. 
		   Useful to fill in information into gaps along a diffusion tensor. */
		template <typename T>
		static BufferedImage<T> inpaintAniso(
			RawImage<T>& data,
			RawImage<T>& weight,
			RawImage<Tensor3x3<float>>& D,
			RawImage<T>& initial,
			double step, int iterations, int num_threads,
			float minVal = 0, float maxVal = 0, bool 
			L1_regional = false,
			bool inflation = false);
		
		/* Much faster version of the above.
		   Works at progressively larger scales indicated in 'scales'.
		   The last value in 'scales' should be 1.0. to obtain an output at the original scale.*/
		template <typename T>
		static BufferedImage<T> inpaintAnisoMultiscale(
			BufferedImage<T>& data,
			BufferedImage<T>& weight,
			BufferedImage<Tensor3x3<float>>& D,
			int num_threads,
			std::vector<double> steps,
			std::vector<double> scales,
			std::vector<int> iterations,
			float minVal = 0, float maxVal = 0, 
			bool L1_regional = false,
			bool inflation = false);
		
		/* Isotropic versions of the above: */
		
		template <typename T>
		static BufferedImage<T> inpaintIso(
			RawImage<T>& data,
			RawImage<T>& weight,
			float D,
			RawImage<T>& initial,
			int iterations, int num_threads);
		
		template <typename T>
		static BufferedImage<T> inpaintIsoMultiscale(
			BufferedImage<T>& data,
			BufferedImage<T>& weight,
			float D,
			int num_threads,
			std::vector<double> scales,
			std::vector<int> iterations);
};

template <typename T>
BufferedImage<T> Diffusion :: diffuse(
		const RawImage<T>& data,
		const RawImage<Tensor3x3<float>>& D,
		float step, int iterations, int num_threads)
{
	const int w = data.xdim;
	const int h = data.ydim;
	const int d = data.zdim;
	
	BufferedImage<gravis::t3Vector<T>> grad(w,h,d);
	BufferedImage<T> u_front_act(w,h,d), u_back_act(w,h,d);
	
	u_back_act = data;
	
	RawImage<T> u_front = u_front_act.getRef();
	RawImage<T> u_back = u_back_act.getRef();
	
	for (int it = 0; it < iterations; it++)
	{
		if (it > 0 && it % 1000 == 0) std::cout << it << std::endl;
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			gravis::f3Vector g;
			
			g.x = x < w-1? u_back(x+1, y, z) - u_back(x,y,z) : u_back(x,y,z) - u_back(x-1,y,z);
			g.y = y < h-1? u_back(x, y+1, z) - u_back(x,y,z) : u_back(x,y,z) - u_back(x,y-1,z);
			g.z = z < d-1? u_back(x, y, z+1) - u_back(x,y,z) : u_back(x,y,z) - u_back(x,y,z-1);
			
			gravis::t3Matrix<T> A = D(x,y,z).toMatrix();
			
			grad(x,y,z) = A * g;
		}
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			T div(0);
			
			if (x > 0) div += grad(x,y,z).x - grad(x-1,y,z).x;
			if (y > 0) div += grad(x,y,z).y - grad(x,y-1,z).y;
			if (z > 0) div += grad(x,y,z).z - grad(x,y,z-1).z;
						
			u_front(x,y,z) = u_back(x,y,z) + step * div;
		}
		
		u_front.swapWith(u_back);
	}
	
	return BufferedImage<T>(u_back);
}


template <typename T>
BufferedImage<T> Diffusion :: diffuseUpward(
		const RawImage<T>& data,
		const RawImage<Tensor3x3<float>>& D,
		float step, int iterations, int num_threads)
{
	const int w = data.xdim;
	const int h = data.ydim;
	const int d = data.zdim;
	
	BufferedImage<gravis::t3Vector<T>> grad(w,h,d);
	BufferedImage<T> u_front_act(w,h,d), u_back_act(w,h,d);
	
	u_back_act = data;
	
	RawImage<T> u_front = u_front_act.getRef();
	RawImage<T> u_back = u_back_act.getRef();
	
	for (int it = 0; it < iterations; it++)
	{
		if (it > 0 && it % 1000 == 0) std::cout << it << std::endl;
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			gravis::f3Vector g;
			
			g.x = x < w-1? u_back(x+1, y, z) - u_back(x,y,z) : u_back(x,y,z) - u_back(x-1,y,z);
			g.y = y < h-1? u_back(x, y+1, z) - u_back(x,y,z) : u_back(x,y,z) - u_back(x,y-1,z);
			g.z = z < d-1? u_back(x, y, z+1) - u_back(x,y,z) : u_back(x,y,z) - u_back(x,y,z-1);
			
			gravis::t3Matrix<T> A = D(x,y,z).toMatrix();
			
			grad(x,y,z) = A * g;
		}
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			T div(0);
			
			if (x > 0) div += grad(x,y,z).x - grad(x-1,y,z).x;
			if (y > 0) div += grad(x,y,z).y - grad(x,y-1,z).y;
			if (z > 0) div += grad(x,y,z).z - grad(x,y,z-1).z;
					
			if (div > 0)
			{
				u_front(x,y,z) = u_back(x,y,z) + step * div;
			}
			else
			{
				u_front(x,y,z) = u_back(x,y,z);
			}
		}
		
		u_front.swapWith(u_back);
	}
	
	return BufferedImage<T>(u_back);
}


template <typename T>
BufferedImage<T> Diffusion :: inpaintAniso(
		RawImage<T>& data,
		RawImage<T>& weight,
		RawImage<Tensor3x3<float>>& D,
		RawImage<T>& initial,
		double step, int iterations, int num_threads,
		float minVal, float maxVal, bool L1_regional, bool inflation)
{
	const int w = data.xdim;
	const int h = data.ydim;
	const int d = data.zdim;
	
	std::cout << "data: " << data.getSizeString() << std::endl;
	std::cout << "weight: " << weight.getSizeString() << std::endl;
	std::cout << "D: " << D.getSizeString() << std::endl;
	std::cout << "initial: " << initial.getSizeString() << std::endl;
	
	BufferedImage<gravis::t3Vector<T>> grad(w,h,d);
	BufferedImage<T> u_front_act(w,h,d), u_back_act(w,h,d);
	
	grad.fill(gravis::t3Vector<T>(0,0,0));
	
	u_back_act = BufferedImage<T>(initial);
	
	RawImage<T> u_front = u_front_act.getRef();
	RawImage<T> u_back = u_back_act.getRef();
	
	for (int it = 0; it < iterations; it++)
	{
		if (it % 100 == 0) std::cout << it << std::endl;
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			gravis::f3Vector g;
			
			g.x = x < w-1? u_back(x+1, y, z) - u_back(x,y,z) : u_back(x-1, y, z) - u_back(x,y,z);
			g.y = y < h-1? u_back(x, y+1, z) - u_back(x,y,z) : u_back(x, y-1, z) - u_back(x,y,z);
			g.z = z < d-1? u_back(x, y, z+1) - u_back(x,y,z) : u_back(x, y, z-1) - u_back(x,y,z);
			
			gravis::t3Matrix<T> A = D(x,y,z).toMatrix();
			
			grad(x,y,z) = A * g;
		}
		
		/*Image<T> divDebug(w,h,d);
		Image<T> gxDebug(w,h,d), gyDebug(w,h,d), gzDebug(w,h,d);*/
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			T div(0);
			
			if (x > 0) div += grad(x,y,z).x - grad(x-1,y,z).x;
			else div += 2 * grad(x,y,z).x;
			
			if (y > 0) div += grad(x,y,z).y - grad(x,y-1,z).y;
			else div += 2 * grad(x,y,z).y;
			
			if (z > 0) div += grad(x,y,z).z - grad(x,y,z-1).z;
			else div += 2 * grad(x,y,z).z;
						
			const T u0 = u_back(x,y,z);
			const T a = weight(x,y,z);
			
			T e;
			
			if (L1_regional)
			{
				e = data(x,y,z) * u0 < 0.0? -u0 / std::abs(u0) : 0.0;
			}
			else
			{
				e = data(x,y,z) - u0;
			}
			
			
			T u1 = u0 + step * (div + a * e);
			
			if (inflation)
			{
				const double inflation_rate = 1.0;
				u1 += step * inflation_rate * u0;
			}
			
			u_front(x,y,z) = u1;
			
			if (minVal != maxVal)
			{
				if (u_front(x,y,z) > maxVal)
				{
					u_front(x,y,z) = maxVal;
				}
				
				if (u_front(x,y,z) < minVal)
				{
					u_front(x,y,z) = minVal;
				}
			}
		}
		
		/*if (it < 500 || (it < 1000 && it % 100 == 0) || (it % 1000 == 0)) 
		{
			divDebug.write("debug/it"+ZIO::itoa(it)+"_divDebug.mrc");
			u_front.write("debug/it"+ZIO::itoa(it)+"_u.mrc");
			data.write("debug/it"+ZIO::itoa(it)+"_data.mrc");
			
			gxDebug.write("debug/it"+ZIO::itoa(it)+"_dgx.mrc");
			gyDebug.write("debug/it"+ZIO::itoa(it)+"_dgy.mrc");
			gzDebug.write("debug/it"+ZIO::itoa(it)+"_dgz.mrc");
		}*/
		
		u_front.swapWith(u_back);
	}
	
	return BufferedImage<T>(u_back);
}

template <typename T>
BufferedImage<T> Diffusion :: inpaintAnisoMultiscale(
		BufferedImage<T>& data,
		BufferedImage<T>& weight,
		BufferedImage<Tensor3x3<float>>& D,
		int num_threads,
		std::vector<double> steps,
		std::vector<double> scales,
		std::vector<int> iterations,
		float minVal, float maxVal, bool L1_regional, bool inflation)
{
	const int w = data.xdim;
	const int h = data.ydim;
	const int d = data.zdim;
	const int sc = steps.size();
			
	BufferedImage<T> initial;
	
	for (int s = 0; s < sc-1; s++)
	{
		const int w1 = w / scales[s];
		const int h1 = h / scales[s];
		const int d1 = d / scales[s];
		
		BufferedImage<T> dataScaled, weightScaled;
		BufferedImage<Tensor3x3<T>> DScaled;
		
		dataScaled = Resampling::downsampleGaussFilt_3D_full(data, w1, h1, d1, scales[s]);
		weightScaled = Resampling::downsampleGaussFilt_3D_full(weight, w1, h1, d1, scales[s]);
		DScaled = StructureTensor::downsample(D, scales[s]);
		
		if (s == 0) initial = dataScaled;
		
		initial = inpaintAniso(
					dataScaled, weightScaled, DScaled, initial, 
					steps[s], 
					iterations[s], 
					w1 > 3*num_threads? num_threads : 1,
					minVal, maxVal, L1_regional, inflation);
		
		initial.write("debug/filled_scale-"+ZIO::itoa(s)+".mrc");
		
		const int w2 = w / scales[s+1];
		const int h2 = h / scales[s+1];
		const int d2 = d / scales[s+1];
		
		initial = Resampling::upsampleLinear_3D_full(
					initial, scales[s] / scales[s+1], w2, h2, d2);
	}
	
	return inpaintAniso(
			data, weight, D, initial, 
			steps[sc-1], 
			iterations[sc-1], 
			data.xdim > 3*num_threads? num_threads : 1,
			minVal, maxVal, L1_regional, inflation);
}

template <typename T>
BufferedImage<T> Diffusion :: inpaintIso(
		RawImage<T>& data,
		RawImage<T>& weight,
		float D,
		RawImage<T>& initial,
		int iterations, int num_threads)
{
	const int w = data.xdim;
	const int h = data.ydim;
	const int d = data.zdim;
	
	BufferedImage<gravis::t3Vector<T>> grad(w,h,d);
	BufferedImage<T> u_front_act(w,h,d), u_back_act(w,h,d);
	
	grad.fill(gravis::t3Vector<T>(0,0,0));
	
	u_back_act = BufferedImage<T>(initial);
	
	RawImage<T> u_front = u_front_act.getRef();
	RawImage<T> u_back = u_back_act.getRef();
	
	for (int it = 0; it < iterations; it++)
	{
		if (it % 100 == 0) std::cout << it << std::endl;
		
		#pragma omp parallel for num_threads(num_threads)
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for	(int x = 0; x < w; x++)
		{
			const float u0 = u_back(x,y,z);
			const float wg = weight(x,y,z);
			
			const float mu = std::max(1.f - wg, 0.f);
			
			float num = 0.f, denom = 0.f;
			
			if (x > 0)
			{
				num += u_back(x-1,y,z) - u0;
				denom += 1.f;
			}
			
			if (x < w-1)
			{
				num += u_back(x+1,y,z) - u0;
				denom += 1.f;
			}
			
			if (y > 0)
			{
				num += u_back(x,y-1,z) - u0;
				denom += 1.f;
			}
			
			if (y < h-1)
			{
				num += u_back(x,y+1,z) - u0;
				denom += 1.f;
			}
			
			if (z > 0)
			{
				num += u_back(x,y,z-1) - u0;
				denom += 1.f;
			}
			
			if (z < d-1)
			{
				num += u_back(x,y,z+1) - u0;
				denom += 1.f;
			}
						
			num  =  mu * num   + wg * (data(x,y,z) - u0);
			denom = mu * denom + wg;
			
			u_front(x,y,z) = u0 + num / denom;
		}
		
		u_front.swapWith(u_back);
	}
	
	return BufferedImage<T>(u_back);
}

template <typename T>
BufferedImage<T> Diffusion :: inpaintIsoMultiscale(
		BufferedImage<T>& data,
		BufferedImage<T>& weight,
		float D,
		int num_threads,
		std::vector<double> scales,
		std::vector<int> iterations)
{
	const int w = data.xdim;
	const int h = data.ydim;
	const int d = data.zdim;
	const int sc = scales.size();
			
	BufferedImage<T> initial;
	
	for (int s = 0; s < sc-1; s++)
	{
		const int w1 = w / scales[s];
		const int h1 = h / scales[s];
		const int d1 = d / scales[s];
		
		BufferedImage<T> dataScaled, weightScaled;
		
		T volScale = scales[s] * scales[s] * scales[s];
		
		dataScaled = Resampling::downsampleGaussFilt_3D_full(data, w1, h1, d1, scales[s]);
		weightScaled = volScale * Resampling::downsampleGaussFilt_3D_full(weight, w1, h1, d1, scales[s]);
		
		if (s == 0) initial = dataScaled;
		
		initial = inpaintIso(
					dataScaled, weightScaled, D, initial, 
					iterations[s], 
					w1 > 3*num_threads? num_threads : 1);
		
		initial.write("debug/filled_scale-"+ZIO::itoa(s)+".mrc");
				
		const int w2 = w / scales[s+1];
		const int h2 = h / scales[s+1];
		const int d2 = d / scales[s+1];
		
		initial = Resampling::upsampleLinear_3D_full(
					initial, scales[s] / scales[s+1], w2, h2, d2);
	}
	
	return inpaintIso(
			data, weight, D, initial, 
			iterations[sc-1], 
			data.xdim > 3*num_threads? num_threads : 1);
}

#endif
