#include "membrane_segmentation.h"
#include <src/jaz/segmentation/diffusion.h>
#include <src/jaz/segmentation/diffusion_tensors.h>

using namespace gravis;


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
