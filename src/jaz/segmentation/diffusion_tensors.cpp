#include "diffusion_tensors.h"

using namespace gravis;


BufferedImage<Tensor3x3<float> > DiffusionTensors::plateDiffusion(
		const RawImage<Tensor3x3<float>>& J, double sigma, double alpha)
{
	BufferedImage<Tensor3x3<float>> D(J.xdim, J.ydim, J.zdim);
	
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	const double sigma2 = sigma * sigma;
		
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		d3Matrix Q;
		d3Vector Lv;
		
		Tensor3x3<double> Jd(J(x,y,z).toMatrix());
		
		Jd.diagonalize(Lv, Q);
		
		float tau0 = alpha + (1-alpha) * exp(-0.5 * Lv[0] * Lv[0] / sigma2);
		
		d3Matrix L(
			tau0, 0, 0, 
			0, 1, 0,
			0, 0, 1);
		
		d3Matrix Qt = Q;
		Qt.transpose();
		
		d3Matrix Dm = Q * L * Qt;
		
		D(x,y,z) = Tensor3x3<float>(Dm);
	}
	
	return D;
}

BufferedImage<Tensor3x3<float> > DiffusionTensors::weightedPlateDiffusion(
		const RawImage<Tensor3x3<float>>& J, 
		const RawImage<float>& edges,
		double sigma, double gamma, double alpha)
{
	BufferedImage<Tensor3x3<float>> D(J.xdim, J.ydim, J.zdim);
	
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	const double sigma2 = sigma * sigma;
		
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		d3Matrix Q;
		d3Vector Lv;
		
		Tensor3x3<double> Jd(J(x,y,z).toMatrix());
		
		Jd.diagonalize(Lv, Q);
		
		float e = edges(x,y,z);
		
		float tau0 = alpha + (1-alpha) * exp(-0.5 * Lv[0] * Lv[0] / sigma2);
		float tau = e * tau0 + (1 - e) * gamma;
		
		d3Matrix L(
			tau, 0, 0, 
			0, 1, 0,
			0, 0, 1);
		
		d3Matrix Qt = Q;
		Qt.transpose();
		
		d3Matrix Dm = Q * L * Qt;
		
		D(x,y,z) = Tensor3x3<float>(Dm);
	}
	
	return D;
}

BufferedImage<Tensor3x3<float> > DiffusionTensors::membraneDiffusion(
		const RawImage<Tensor3x3<float>>& J, 
		const RawImage<float>& edges)
{
	BufferedImage<Tensor3x3<float>> D(J.xdim, J.ydim, J.zdim);
	
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
			
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		d3Matrix Q;
		d3Vector Lv;
		
		Tensor3x3<double> Jd(J(x,y,z).toMatrix());
		
		Jd.diagonalize(Lv, Q);
		
		float e = edges(x,y,z);
		
		d3Matrix L(
			0, 0, 0, 
			0, e, 0,
			0, 0, e);
		
		d3Matrix Qt = Q;
		Qt.transpose();
		
		d3Matrix Dm = Q * L * Qt;
		
		D(x,y,z) = Tensor3x3<float>(Dm);
	}
	
	return D;
}

BufferedImage<Tensor3x3<float> > DiffusionTensors::antiPlateDiffusion(
		const RawImage<Tensor3x3<float>>& J, double alpha)
{
	BufferedImage<Tensor3x3<float>> D(J.xdim, J.ydim, J.zdim);
	
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
		
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		d3Matrix Q;
		d3Vector Lv;
		
		Tensor3x3<double> Jd(J(x,y,z).toMatrix());
		
		Jd.diagonalize(Lv, Q);
		
		d3Matrix L(
			1, 0, 0, 
			0, alpha, 0,
			0, 0, alpha);
		
		d3Matrix Qt = Q;
		Qt.transpose();
		
		d3Matrix Dm = Q * L * Qt;
		
		D(x,y,z) = Tensor3x3<float>(Dm);
	}
	
	return D;
}
