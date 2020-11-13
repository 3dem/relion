#include "membrane_segmentation.h"
#include <src/jaz/segmentation/diffusion.h>
#include <src/jaz/segmentation/diffusion_tensors.h>

using namespace gravis;


BufferedImage<float> MembraneSegmentation::constructMembraneKernel(
		int w, int h, int d,
		double falloff, double kernel_width, double spacing, double ratio, double depth,
		double angle)
{
	BufferedImage<float> kernel(w,h,d);
	kernel.fill(0.f);
	
	const double s2x = 2 * kernel_width * kernel_width;
	const double s2y = 2 * falloff * falloff;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		//const double xx = x < w/2? x : x - w;
		//const double yy = (y < h/2? y : y - h) - depth * spacing;

		const double x0 = x < w/2? x : x - w;
		const double y0 = y < h/2? y : y - h;

		const double xx = cos(angle) * x0 - sin(angle) * y0;
		const double yy = sin(angle) * x0 + cos(angle) * y0 - depth * spacing;
		
		const double yt = yy / spacing;
		
		double wave;
		
		if (yt < 0.5 && yt > -1.5)
		{
			wave = sin(2*PI*yt);
			
			if (yt < 0)
			{
				wave *= (yt + 3) / 3;
			}
		}
		else
		{
			wave = 0.0;
		}
		
		double decay = yy >= 0 ? exp(-yy*yy/s2y) : -exp(-yy*yy/s2y);
		
		kernel(x,y) = exp(-xx*xx/s2x) * (ratio * wave + decay);
	}
	
	return kernel;
}

BufferedImage<float> MembraneSegmentation::correlateWithMembrane(
		BufferedImage<float>& map,
		double falloff, double kernel_width, double spacing, double ratio, double depth,
		double angle)
{
	BufferedImage<float> kernel = MembraneSegmentation::constructMembraneKernel(
		map.xdim, map.ydim, map.zdim, falloff, kernel_width, spacing, ratio, depth, angle);

	BufferedImage<fComplex> map_FS, kernel_FS, correlation_FS;

	FFT::FourierTransform(map, map_FS);
	FFT::FourierTransform(kernel, kernel_FS);

	correlation_FS.resize(map_FS.xdim, map_FS.ydim, map_FS.zdim);

	for (int z = 0; z < map_FS.zdim; z++)
	for (int y = 0; y < map_FS.ydim; y++)
	for (int x = 0; x < map_FS.xdim; x++)
	{
		correlation_FS(x,y,z) = map_FS(x,y,z) * kernel_FS(x,y,z).conj();
	}

	BufferedImage<float> correlation;

	FFT::inverseFourierTransform(correlation_FS, correlation);

	return correlation;
}

BufferedImage<float> MembraneSegmentation::correlateWithMembraneMultiAngle(
		BufferedImage<float>& map,
		double falloff, double kernel_width, double spacing, double ratio, double depth, double max_tilt, int tilt_steps)
{
	BufferedImage<float> correlation = correlateWithMembrane(
		map, falloff, kernel_width, spacing, ratio, depth);

	for (int t = -tilt_steps/2; t <= tilt_steps/2; t++)
	{
		if (t == 0) continue;

		BufferedImage<float> correlation_t = correlateWithMembrane(
			map, falloff, kernel_width, spacing, ratio, depth,
			2 * t * max_tilt / tilt_steps);

		for (int i = 0; i < correlation.getSize(); i++)
		{
			if (correlation_t[i] > correlation[i])
			{
				correlation[i] = correlation_t[i];
			}
		}
	}

	return correlation;
}

BufferedImage<float> MembraneSegmentation::determineMembraniness(
		const RawImage<float>& tomo,
		RawImage<Tensor3x3<float>>& J,
		double sigma_diff, 
		double lambda_fin,
		double thresh_edge,
		int num_threads,
		int itersAcross,
		int itersAlong,
		float stepSize)
{
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	BufferedImage<Tensor3x3<float>> D = DiffusionTensors::antiPlateDiffusion(J);
			
	BufferedImage<float> diffusedAcross = Diffusion::diffuse(
				tomo, D, stepSize, itersAcross, num_threads);
	
	BufferedImage<float> memb = diffusedAcross - tomo;
	
	D = DiffusionTensors::plateDiffusion(J, sigma_diff);
	
	BufferedImage<float> diffusedAlong = Diffusion::diffuse(
				memb, D, stepSize, itersAlong, num_threads);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		const float m0 = memb(x,y,z);
		const float m1 = diffusedAlong(x,y,z);
		
		memb(x,y,z) = m1 > 0? m0 * m1 : 0;
	}
	
	memb = Diffusion::diffuse(memb, D, stepSize, itersAcross, num_threads);
	
	/*Image<float> maxAcross = Diffusion::diffuseUpward(
				memb, D, stepSize, itersAcross, num_threads);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		const float f = memb(x,y,z);
		const float m = maxAcross(x,y,z);
		
		memb(x,y,z) = (f > 0 && m > 0)? f / (m + lambda_fin) : 0;
		
		if (memb(x,y,z) < thresh_edge)
		{
			memb(x,y,z) = 0.0;
		}
		else
		{
			memb(x,y,z) = (memb(x,y,z) - thresh_edge) / (1 - thresh_edge);
		}
	}*/
	
	return memb;
}

BufferedImage<float> MembraneSegmentation::findSilhouette3D(
		const RawImage<float>& membrane, 
		RawImage<Tensor3x3<float>>& J, 
		d4Matrix view, 
		double kappa)
{
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	BufferedImage<float> out(w,h,d);
	
	d3Vector dir_z = d3Vector(view(2,0), view(2,1), view(2,2)).normalize();
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		f3Matrix Q;
		f3Vector L;
		
		J(x,y,z).diagonalize(L, Q);
		
		const d3Vector dir_n = d3Vector(Q(0,0),Q(1,0),Q(2,0)).normalize();
		
		const double zn = dir_z.dot(dir_n);
		const double xyn = sqrt(1.0 - zn*zn);
		
		out(x,y,z) = pow(xyn, kappa) * membrane(x,y,z);
	}
	
	return out;
}

BufferedImage<float> MembraneSegmentation::softmaxMembraneDist(
		const RawImage<float>& membrane, 
		float lambda)
{
	const int w = membrane.xdim;
	const int wh = w/2 + 1;
	const int h = membrane.ydim;
	const int d = membrane.zdim;
		
	BufferedImage<float> kernel(w,h,d);
		
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		const double xx = x < w/2? x : x - w;
		const double yy = y < h/2? y : y - h;
		const double zz = z < d/2? z : z - d;
		
		const double r = sqrt(xx*xx + yy*yy + zz*zz);
		
		kernel(x,y,z) = exp(-r/lambda);
	}
		
	BufferedImage<float> membraneCp = membrane;
	
	BufferedImage<fComplex> kernelFS, membraneFS;
				
	FFT::FourierTransform(kernel, kernelFS, FFT::Both);
	FFT::FourierTransform(membraneCp, membraneFS, FFT::Both);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < wh; x++)
	{
		membraneFS(x,y,z) *= kernelFS(x,y,z);
	}
	
	FFT::inverseFourierTransform(membraneFS, membraneCp, FFT::Both);
	
	return membraneCp;	
}
