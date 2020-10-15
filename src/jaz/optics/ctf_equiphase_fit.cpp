#include "ctf_equiphase_fit.h"

using namespace gravis;

CtfEquiphaseFit::CtfEquiphaseFit(
		const BufferedImage<float> &spectrum, 
		double pixelSize, 
		double voltage, 
		double Cs, 
		double Q0, 
		double mapScale,
		int num_threads, 
		int k0, int k1)
	:
	spectrum(spectrum),
	num_threads(num_threads), k0(k0), k1(k1),
	pixelSize(pixelSize), voltage(voltage), Cs(Cs), Q0(Q0),
	mapScale(mapScale)
{
	exampleCtf.setValues(0, 0, 0, voltage, Cs, Q0, 0.0, 1.0, 0.0);
}

void CtfEquiphaseFit::averageAlongIso(
		const CTF& ctf,
		std::vector<double>& accum, 
		std::vector<double>& weight,
		int boundary) const
{
	const int rad = accum.size();
	const int s = spectrum.ydim;
	const int sh = spectrum.xdim;
	
	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;
	
	for (int r = 0; r < rad; r++)
	{
		accum[r] = 0.0;
		weight[r] = 0.0;
	}
	
	const double p2a = 1.0 / (double)(s * pixelSize);
	
	for (int yi = b0; yi < s-b1;  yi++)
	for (int xi = b0; xi < sh; xi++)
	{
		const double kx = xi;
		const double ky = yi < s/2? yi : yi - s;
		
		const double r2 = -mapScale * ctf.getLowOrderGamma(p2a * kx, p2a * ky);
		
		const double r = r2 > 0.0? sqrt(r2) : 0.0;
		
		const int r0 = (int) r;
		const int r1 = r0 + 1;
		const double dr = r - r0;
		
		const float val = spectrum(xi,yi);
		
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

double CtfEquiphaseFit::compareWithExpansion(
		const CTF& ctf,
		const std::vector<double>& accum,
		int boundary) const
{
	const int rad = accum.size();
	const int s = spectrum.ydim;
	const int sh = spectrum.xdim;
	
	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;
	
	const double p2a = 1.0 / (double)(s * pixelSize);
	
	double out(0.0);
	
	for (int yi = b0; yi < s-b1;  yi++)
	for (int xi = b0; xi < sh; xi++)
	{
		const double kx = xi;
		const double ky = yi < s/2? yi : yi - s;
		
		const double r2f = kx * kx + ky * ky;
		
		if (r2f < k0*k0 || r2f > k1*k1) continue;
		
		const double r2 = -mapScale * ctf.getLowOrderGamma(p2a * kx, p2a * ky);
		
		const double r = r2 > 0.0? sqrt(r2) : 0.0;
		
		const int r0 = r < rad - 2? (int) r : rad - 2;
		const int r1 = r0 + 1;
		const double dr = r < rad - 2? r - r0 : 1.0;
		
		const float v = spectrum(xi,yi);
		const double u = (1.0 - dr) * accum[r0] + dr * accum[r1];
					
		out += (u - v) * (u - v);
	}
	
	return out;
}

BufferedImage<double> CtfEquiphaseFit::computeExpansion(
		const CTF& ctf,
		const std::vector<double>& accum,
		int boundary) const
{
	const int rad = accum.size();
	const int s = spectrum.ydim;
	const int sh = spectrum.xdim;
	
	const int b0 = boundary;
	const int b1 = boundary < 2? 0 : boundary - 1;
	
	const double p2a = 1.0 / (double)(s * pixelSize);
	
	BufferedImage<double> out(sh,s);
	
	for (int yi = b0; yi < s-b1;  yi++)
	for (int xi = b0; xi < sh; xi++)
	{
		const double kx = xi;
		const double ky = yi < s/2? yi : yi - s;
		
		if (kx*kx + ky*ky < k0*k0) continue;
		
		const double r2 = -mapScale * ctf.getLowOrderGamma(p2a * kx, p2a * ky);
		
		const double r = r2 > 0.0? sqrt(r2) : 0.0;
		
		const int r0 = r < rad - 2? (int) r : rad - 2;
		const int r1 = r0 + 1;
		const double dr = r < rad - 2? r - r0 : 1.0;
		
		out(xi,yi) = (1.0 - dr) * accum[r0] + dr * accum[r1];
	}
	
	return out;
}

double CtfEquiphaseFit::f(const std::vector<double>& x, void* tempStorage) const
{
	const int sh = spectrum.xdim;
	
	const int rad = 2 * sh;
		
	CTF ctf = paramsToCtf(x);
		
	std::vector<double> accum(rad), weight(rad);
	
	averageAlongIso(ctf, accum, weight, 2);
	
	return compareWithExpansion(ctf, accum, 2);
}

void CtfEquiphaseFit::report(int iteration, double cost, const std::vector<double> &x) const
{
	std::cout << iteration << ": " << cost << " @  " 
			  << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
}

CTF CtfEquiphaseFit::paramsToCtf(const std::vector<double> &x) const
{
	CTF ctf = exampleCtf;
	ctf.setDefocusMatrix(x[0], x[1], x[2]);
	return ctf;
}
