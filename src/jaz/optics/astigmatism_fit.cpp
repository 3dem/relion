#include "astigmatism_fit.h"

AstigmatismFit::AstigmatismFit(
		const std::vector<BufferedImage<float>>& spectra, 
		int num_threads, 
		int k0,
		int k1,
		double Cs_px )
	:	
	spectra(spectra),
	num_threads(num_threads),
	k0(k0),
	k1(k1),
	Cs_px(Cs_px)
{		
	std::cout << spectra.size() << " frames." << std::endl;
}

void AstigmatismFit::averageAlongIso(
		int f, 
		double alpha, double beta, 
		std::vector<double>& accum, 
		std::vector<double>& weight,
		int boundary) const
{
	const int rad = accum.size();
	const int s = spectra[0].ydim;
	const int sh = spectra[0].xdim;
	
	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;
	
	for (int r = 0; r < rad; r++)
	{
		accum[r] = 0.0;
		weight[r] = 0.0;
	}
	
	for (int yi = b0; yi < s-b1;  yi++)
	for (int xi = b0; xi < sh; xi++)
	{
		const double kx = xi;
		const double ky = yi < s/2? yi : yi - s;
		
		const double r2 = 
				kx * kx + ky * ky 
				+ alpha * (kx * kx - ky * ky)
				+ 2.0 * beta * kx * ky;
		
		const double r = r2 > 0.0? sqrt(r2) : 0.0;
		
		const int r0 = (int) r;
		const int r1 = r0 + 1;
		const double dr = r - r0;
		
		const float val = spectra[f](xi,yi);
		
		if (r0 < rad)
		{
			accum[r0] += (1.0 - dr) * val;
			weight[r0] += (1.0 - dr);
		}
		
		if (r1 < rad)
		{
			accum[r1] += dr * val;
			weight[r1] += dr;
		}
	}
	
	for (int r = 0; r < rad; r++)
	{
		if (weight[r] > 0.0)
		{
			accum[r] /= weight[r];
		}
	}
}

double AstigmatismFit::compareWithExpansion(
		int f, 
		double alpha, double beta, 
		const std::vector<double>& accum,
		int boundary) const
{
	const int rad = accum.size();
	const int s = spectra[0].ydim;
	const int sh = spectra[0].xdim;
	
	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;
	
	double out(0.0);
	
	for (int yi = b0; yi < s-b1;  yi++)
	for (int xi = b0; xi < sh; xi++)
	{
		const double kx = xi;
		const double ky = yi < s/2? yi : yi - s;
		
		const double r2f = kx * kx + ky * ky;
		
		if (r2f < k0*k0 || r2f > k1*k1) continue;
		
		const double r2 = 
				r2f
				+ alpha * (kx * kx - ky * ky)
				+ 2.0 * beta * kx * ky;
		
		const double r = r2 > 0.0? sqrt(r2) : 0.0;
		
		const int r0 = r < rad - 2? (int) r : rad - 2;
		const int r1 = r0 + 1;
		const double dr = r < rad - 2? r - r0 : 1.0;
		
		const float v = spectra[f](xi,yi);
		const double u = (1.0 - dr) * accum[r0] + dr * accum[r1];
					
		out += (u - v) * (u - v);
	}
	
	return out;
}

BufferedImage<double> AstigmatismFit::computeExpansion(		
		int f, 
		double alpha, double beta, 
		const std::vector<double>& accum,
		int boundary) const
{
	const int rad = accum.size();
	const int s = spectra[0].ydim;
	const int sh = spectra[0].xdim;
	
	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;
	
	BufferedImage<double> out(sh,s);
	
	for (int yi = b0; yi < s-b1;  yi++)
	for (int xi = b0; xi < sh; xi++)
	{
		const double kx = xi;
		const double ky = yi < s/2? yi : yi - s;
		
		if (kx*kx + ky*ky < k0*k0) continue;
		
		const double r2 = 
				kx * kx + ky * ky 
				+ alpha * (kx * kx - ky * ky)
				+ 2.0 * beta * kx * ky;
		
		const double r = r2 > 0.0? sqrt(r2) : 0.0;
		
		const int r0 = r < rad - 2? (int) r : rad - 2;
		const int r1 = r0 + 1;
		const double dr = r < rad - 2? r - r0 : 1.0;
		
		out(xi,yi) = (1.0 - dr) * accum[r0] + dr * accum[r1];
	}
	
	return out;
}

double AstigmatismFit::f(const std::vector<double>& x, void* tempStorage) const
{
	const int fc = spectra.size();
	const int sh = spectra[0].xdim;
	
	const int rad = 2 * sh;
	
	double out(0.0);
	
	const double alpha = x[0];
	const double beta = x[1];
	
	for (int f = 0; f < fc; f++)
	{
		std::vector<double> accum(rad), weight(rad);
		
		averageAlongIso(f, alpha, beta, accum, weight, 2);
		
		out += compareWithExpansion(f, alpha, beta, accum, 2);
	}
	
	return out;
}

void AstigmatismFit::report(int iteration, double cost, const std::vector<double> &x) const
{
	std::cout << iteration << ": " << cost << " @  " << x[0] << ", " << x[1] << std::endl;
}
