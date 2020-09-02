#include "spectral_ctf_cost.h"
#include <src/jaz/image/interpolation.h>

using namespace gravis;

#define DEFOCUS 0
#define ALPHA 1
#define BETA 2
#define SCALE 3
#define BFAC 4
#define ZETA 5

#define PARAM_COUNT 6

#define Z_SCALE 1e2
#define B_SCALE 1e-4


SpectralCtfCost::SpectralCtfCost()
:
	num_threads(1), k0(0), k1(0),
	pixelSize(0), voltage(0), Cs(0), Q0(0)
{	
}

SpectralCtfCost::SpectralCtfCost(
		RawImage<float> spectrum, 
		RawImage<float> background,
		double pixelSize, 
		double voltage, 
		double Cs, 
		double Q0, 
		int num_threads, 
		int k0, int k1)
	:
	spectrum(spectrum),
	background(background),
	num_threads(num_threads), k0(k0), k1(k1),
	pixelSize(pixelSize), voltage(voltage), Cs(Cs), Q0(Q0)
{	
	const double phase_shift = 0.0;
	
	// taken from Relion's CTF class:
	
	double local_Cs = Cs * 1e7;
    double local_kV = voltage * 1e3;

    lambda = 12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));
	
	const double p2a = 1.0 / (pixelSize * spectrum.ydim);
	const double p2a2 = p2a * p2a;

    R1 =  PI * lambda * p2a2;
    R2 = (PI / 2) * local_Cs * lambda * lambda * lambda * p2a2 * p2a2;
    R3_5 = atan(Q0/sqrt(1-Q0*Q0)) + DEG2RAD(phase_shift);
}

double SpectralCtfCost::f(const std::vector<double>& x, void* tempStorage) const
{
	const double dZ    = Z_SCALE * x[DEFOCUS];
	const double alpha = Z_SCALE * x[ALPHA];
	const double beta  = Z_SCALE * x[BETA];
	const double scale = x[SCALE];
	const double B4px  = B_SCALE * 4.0 * x[BFAC];
	const double zeta  = x[ZETA];
	
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	const double Axx = -dZ + alpha;
	const double Axy = beta;
	const double Ayy = -dZ - alpha;
		
	const double k0sq = k0 * k0;
	const double k1sq = k1 * k1;
	
	const int b0 = 2;
	const int b1 = 1;
	
	
	double out(0.0);
	
	for (int yi = b0; yi < h-b1;  yi++)
	for (int xi = b0; xi < wh; xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
		
		if (k2 < k0sq || k2 > k1sq) continue;
		
		const double k4 = k2 * k2;
		
		double gamma = R1 * (Axx*kx*kx + 2.0*Axy*kx*ky + Ayy*ky*ky) + R2 * k4 - R3_5;
		
		const double ctf = -scale * exp(-B4px * k2) * cos(2.0 * gamma);
		
		const double y = zeta * background(xi,yi) + ctf;
		
		const double e = y - spectrum(xi,yi);
		
		out += e*e;
	}
	
	return out;
}

void SpectralCtfCost::grad(
	const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const
{
	const double dZ    = Z_SCALE * x[DEFOCUS];
	const double alpha = Z_SCALE * x[ALPHA];
	const double beta  = Z_SCALE * x[BETA];
	const double scale = x[SCALE];
	const double B4px  = B_SCALE * 4.0 * x[BFAC];
	const double zeta  = x[ZETA];
	
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	const double Axx = -dZ + alpha;
	const double Axy = beta;
	const double Ayy = -dZ - alpha;
		
	const double k0sq = k0 * k0;
	const double k1sq = k1 * k1;
	
	const int b0 = 2;
	const int b1 = 1;
		
	
	for (int i = 0; i < PARAM_COUNT; i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int yi = b0; yi < h-b1; yi++)
	for (int xi = b0; xi < wh;   xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
		
		if (k2 < k0sq || k2 > k1sq) continue;
		
		const double k4 = k2 * k2;
		
		const double astig_defocus = Axx*kx*kx + 2.0*Axy*kx*ky + Ayy*ky*ky;
		
		double gamma = R1 * astig_defocus + R2 * k4 - R3_5;
		
		const double envelope = exp(-B4px * k2);
		const double wave = -cos(2.0 * gamma);
		
		const double y = zeta * background(xi,yi) + scale * envelope * wave;
		
		const double e = y - spectrum(xi,yi);
		
		// C += e*e;
		
		const double dC_dy = 2.0 * e;
		const double dwave_dgamma = 2.0 * sin(2.0 * gamma);
		const double dC_dgamma = dC_dy * scale * envelope * dwave_dgamma;
		
		const double dgamma_ddZ = -R1 * k2;
		const double dgamma_dalpha = R1 * (kx*kx - ky*ky);
		const double dgamma_dbeta = 2.0 * R1 * kx*ky;
		
		
		gradDest[DEFOCUS] += dC_dgamma * dgamma_ddZ    * Z_SCALE;
		gradDest[ALPHA]   += dC_dgamma * dgamma_dalpha * Z_SCALE;
		gradDest[BETA]    += dC_dgamma * dgamma_dbeta  * Z_SCALE;
		
		gradDest[SCALE]   += dC_dy * envelope * wave;
		gradDest[BFAC]    += dC_dy * scale * wave * envelope * (-k2) * 4.0 * B_SCALE;
		gradDest[ZETA]    += dC_dy * background(xi,yi);
	}
}

void SpectralCtfCost::report(int iteration, double cost, const std::vector<double> &x) const
{
	std::cout << iteration << ": " << cost << " @  [" 
			  << x[SCALE] << ", " << B_SCALE * x[BFAC] << "] * [" 
			  << Z_SCALE * x[DEFOCUS] << ", " 
			  << Z_SCALE * x[ALPHA] << ", " 
			  << Z_SCALE * x[BETA] << "] + [" 
			  << x[ZETA] << "]"
			  << std::endl;
}

std::vector<double> SpectralCtfCost::getInitialParams(double defocus)
{
	std::vector<double>	out(PARAM_COUNT);
	
	out[DEFOCUS] = defocus / Z_SCALE;
	out[ALPHA] = 0.0;
	out[BETA]  = 0.0;
	out[SCALE] = 0.01;
	out[BFAC]  = 0.0;
	out[ZETA]  = 1.0;
	
	return out;
}

BufferedImage<float> SpectralCtfCost::render(const std::vector<double> &x)
{
	const double dZ    =     Z_SCALE * x[DEFOCUS];
	const double alpha =     Z_SCALE * x[ALPHA];
	const double beta  =     Z_SCALE * x[BETA];
	const double scale =               x[SCALE];
	const double B4px  = 4 * B_SCALE * x[BFAC];
	const double zeta  =               x[ZETA];
	
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	const double Axx = -dZ + alpha;
	const double Axy = beta;
	const double Ayy = -dZ - alpha;
	
	
	BufferedImage<float> out(wh,h);
	
	for (int yi = 0; yi < h;  yi++)
	for (int xi = 0; xi < wh; xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
		const double k4 = k2 * k2;
		
		double gamma = R1 * (Axx*kx*kx + 2.0*Axy*kx*ky + Ayy*ky*ky) + R2 * k4 - R3_5;
		const double ctf = -scale * exp(-B4px*k2) * cos(2.0 * gamma);
		
		out(xi,yi) = zeta * background(xi,yi) + ctf;
	}
	
	return out;
}

void SpectralCtfCost::addAlignedSpectrum(
	double z0, const std::vector<double> &x, BufferedImage<double> &accum, BufferedImage<double> &weight)
{
	const double z0_A  =     Z_SCALE * z0;
	const double z1_A  =     Z_SCALE * x[DEFOCUS];
	const double alpha =     Z_SCALE * x[ALPHA];
	const double beta  =     Z_SCALE * x[BETA];
	
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	const double Bxx = alpha;
	const double Bxy = beta;
	const double Byy = -alpha;
	
	
	for (int yi = 0; yi < h;  yi++)
	for (int xi = 0; xi < wh; xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
		const double k4 = k2 * k2;
		
		const double ktBk = Bxx*kx*kx + 2.0*Bxy*kx*ky + Byy*ky*ky;
		
		const double d = R1 * R1 * (ktBk - z0_A * k2) * (ktBk - z0_A * k2)
			- 4.0 * R2 * k4 * (R1 * (z1_A * k2 - ktBk) - R2 * k4);
		
		if (d > 0.0)
		{
			const double cc = (R1 * (z0_A * k2 - ktBk) - sqrt(d)) / (2.0 * R2 * k4);
			
			if (cc > 0.0)
			{
				const double c = sqrt(cc);
				
				const double ckx = c * kx;
				const double cky = c * ky;
				
				Interpolation::insertLinearXY_FftwHalf(
					(double) spectrum(xi,yi), 1.0, ckx, cky, accum, weight);
			}
		}
	}	
}

void SpectralCtfCost::addAlignedRadialAverage(
		double z0, const std::vector<double> &x, 
		std::vector<double> &accum, std::vector<double> &weight)
{
	const double z0_A  =     Z_SCALE * z0;
	const double z1_A  =     Z_SCALE * x[DEFOCUS];
	const double alpha =     Z_SCALE * x[ALPHA];
	const double beta  =     Z_SCALE * x[BETA];
	
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	const double Bxx = alpha;
	const double Bxy = beta;
	const double Byy = -alpha;
	
	
	for (int yi = 0; yi < h;  yi++)
	for (int xi = 0; xi < wh; xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
		const double k4 = k2 * k2;
		
		const double ktBk = Bxx*kx*kx + 2.0*Bxy*kx*ky + Byy*ky*ky;
		
		const double d = R1 * R1 * (ktBk - z0_A * k2) * (ktBk - z0_A * k2)
			- 4.0 * R2 * k4 * (R1 * (z1_A * k2 - ktBk) - R2 * k4);
		
		if (d > 0.0)
		{
			const double cc = (R1 * (z0_A * k2 - ktBk) - sqrt(d)) / (2.0 * R2 * k4);
			
			if (cc > 0.0)
			{
				const double c = sqrt(cc);
				
				const double ckx = c * kx;
				const double cky = c * ky;
				
				const double r = sqrt(ckx * ckx + cky * cky);
				const int r0 = (int) r;
				const int r1 = r0 + 1;
				const double dr = r - r0;
				
				if (r0 < wh)
				{
					accum[r0] += (1.0 - dr) * spectrum(xi,yi);
					weight[r0] += (1.0 - dr);
				}
				
				if (r1 < wh)
				{
					accum[r1] += dr * spectrum(xi,yi);
					weight[r1] += dr;
				}
			}
		}
	}
}


void SpectralCtfCost::addRadialAverage(
		std::vector<double> &accum, std::vector<double> &weight)
{
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	for (int yi = 0; yi < h;  yi++)
	for (int xi = 0; xi < wh; xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
		
		const double r = sqrt(k2);
		const int r0 = (int) r;
		const int r1 = r0 + 1;
		const double dr = r - r0;
		
		if (r0 < wh)
		{
			accum[r0] += (1.0 - dr) * spectrum(xi,yi);
			weight[r0] += (1.0 - dr);
		}
		
		if (r1 < wh)
		{
			accum[r1] += dr * spectrum(xi,yi);
			weight[r1] += dr;
		}
	}
}

void SpectralCtfCost::addRadialAverageFit(
		const std::vector<double> &x, 
		std::vector<double> &accum, 
		std::vector<double> &weight)
{
	const double dZ    = Z_SCALE * x[DEFOCUS];
	const double alpha = Z_SCALE * x[ALPHA];
	const double beta  = Z_SCALE * x[BETA];
	const double scale = x[SCALE];
	const double B4px  = B_SCALE * 4.0 * x[BFAC];
	const double zeta  = x[ZETA];
	
	const int wh = spectrum.xdim;
	const int h  = spectrum.ydim;
	
	const double Axx = -dZ + alpha;
	const double Axy = beta;
	const double Ayy = -dZ - alpha;
		
	const double k0sq = k0 * k0;
	const double k1sq = k1 * k1;
	
	const int b0 = 2;
	const int b1 = 1;
	
	
	for (int yi = b0; yi < h-b1;  yi++)
	for (int xi = b0; xi < wh; xi++)
	{
		const double kx = xi;
		const double ky = yi < h/2? yi : yi - h;
		const double k2 = kx * kx + ky * ky;
				
		const double k4 = k2 * k2;
		
		double gamma = R1 * (Axx*kx*kx + 2.0*Axy*kx*ky + Ayy*ky*ky) + R2 * k4 - R3_5;
		
		const double ctf = -scale * exp(-B4px * k2) * cos(2.0 * gamma);
		
		const double y = zeta * background(xi,yi) + ctf;
		
		const double r = sqrt(k2);
		const int r0 = (int) r;
		const int r1 = r0 + 1;
		const double dr = r - r0;
		
		if (r0 < wh)
		{
			accum[r0] += (1.0 - dr) * y;
			weight[r0] += (1.0 - dr);
		}
		
		if (r1 < wh)
		{
			accum[r1] += dr * y;
			weight[r1] += dr;
		}
	}
}

double SpectralCtfCost::offsetDefocusParam(double x0, double deltaZ_Ang) const
{
	return x0 + deltaZ_Ang / Z_SCALE;
}
