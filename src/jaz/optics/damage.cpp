#include "damage.h"
#include <src/jaz/image/power_spectrum.h>

using namespace gravis;


void Damage::applyWeight(
		RawImage<float>& stack,
		double pixelSize,
		const std::vector<double>& doses,
		int num_threads)
{
	const int w  = stack.xdim;
	const int h  = stack.ydim;
	const int fc = stack.zdim;

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<float> frame = stack.getConstSliceRef(f);
		BufferedImage<fComplex> frame_FS;

		FFT::FourierTransform(frame, frame_FS);

		const int wh = frame_FS.xdim;

		for (int y = 0; y < h;  y++)
		for (int x = 0; x < wh; x++)
		{
			const double xx = x / (pixelSize * w);
			const double yy = (y < h/2? y : y - h) / (pixelSize * h);

			const double k = sqrt(xx*xx + yy*yy);

			frame_FS(x,y) *= getWeight(doses[f], k);
		}

		FFT::inverseFourierTransform(frame_FS, frame);

		stack.getSliceRef(f).copyFrom(frame);
	}
}

void Damage::applyWeight(
		RawImage<fComplex>& stack,
		double pixelSize,
		const std::vector<double>& doses,
		int num_threads)
{
	const int wh = stack.xdim;
	const int h  = stack.ydim;
	const int fc = stack.zdim;
	const int w = 2 * (wh - 1);

	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h;  y++)
		for (int x = 0; x < wh; x++)
		{
			const double xx = x / (pixelSize * w);
			const double yy = (y < h/2? y : y - h) / (pixelSize * h);

			const double k = sqrt(xx*xx + yy*yy);

			stack(x,y,f) *= getWeight(doses[f], k);
		}
	}
}

std::vector<double> Damage::criticalDamageVector(double pixelSize, int boxSize)
{
	const int sh = boxSize / 2 + 1;
	
	const double a = 0.245;
	const double b = -1.665;
	const double c = 2.81;
	
	std::vector<double> out(sh);
	
	for (int r = 0; r < sh; r++)
	{
		const double k = r / (boxSize * pixelSize);
		out[r] = a * pow(k, b) + c;
	}
	
	return out;
}

std::vector<double> Damage::weightVector(double dose, double pixelSize, int boxSize)
{
	std::vector<double> out = criticalDamageVector(pixelSize, boxSize);
	
	const int sh = boxSize / 2 + 1;
	
	for (int r = 0; r < sh; r++)
	{
		out[r] = exp(-0.5 * dose / out[r]);
	}
	
	return out;
}

BufferedImage<float> Damage::weightImage(double dose, double pixelSize, int boxSize)
{
	const double a = 0.245;
	const double b = -1.665;
	const double c = 2.81;
	
	const int s = boxSize;
	const int sh = s/2 + 1;
	
	
	BufferedImage<float> out(sh,s);
	
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		double xx = x;
		double yy = y < s/2? y : y - s;
		double r = sqrt(xx*xx + yy*yy);
		
		const double k = r / (boxSize * pixelSize);
		const double d0 = a * pow(k, b) + c;

		out(x,y) = (float) exp(-0.5 * dose / d0);
	}
	
	return out;
}

BufferedImage<float> Damage::weightStack_GG(
		const std::vector<double> &doses, double pixelSize, int boxSize)
{
	const int fc = doses.size();
	const int s = boxSize;
	const int sh = s/2 + 1;
	
	BufferedImage<float> out(sh,s,fc);
			
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<float> df = weightImage(doses[f], pixelSize, boxSize);
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			out(x,y,f) = df(x,y);
		}
	}
	
	return out;
}

double Damage::getWeight(double dose, double k)
{
	const double a = 0.245;
	const double b = -1.665;
	const double c = 2.81;
	
	return exp(-0.5 * dose / (a * pow(k, b) + c));
}

std::vector<double> Damage::estimateDoseWeights(
	const BufferedImage<float>& tiltSeries, 
	double pixelSize, int tileSize,
	int firstFrame, double k0_Ang, bool diag)
{
	const int fc = tiltSeries.zdim;
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	const int r0dmg_px = s * pixelSize / k0_Ang;
	
	std::vector<std::vector<double>> powSum(fc);
	
	for (int f = 0; f < fc; f++)
	{
		powSum[f] = PowerSpectrum::periodogramAverage1D(tiltSeries, tileSize, 2.0, f, true);
	}
	
	int refFrame(0);
	
	if (firstFrame < 0)
	{
		double maxFramePower(0.0);
		int bestFrame(0);
		
		for (int f = 0; f < fc; f++)
		{
			double framePower(0.0);
			
			for (int r = r0dmg_px; r < sh; r++)
			{
				framePower += powSum[f][r];
			}
			
			if (framePower > maxFramePower)
			{
				maxFramePower = framePower;
				bestFrame = f;
			}
		}
		
		refFrame = bestFrame;
	}
	else
	{
		refFrame = firstFrame;
	}
	
	std::cout << "    estimating dose weighting using frame " << refFrame 
			  << " as reference..." << std::endl;
	
	std::vector<double> Nk = Damage::criticalDamageVector(pixelSize, s);
		
	std::vector<double> frameDose(fc);
	
	std::cout << "    cumulative doses detected per frame [e/Å²]:";
		
	for (int f = 0; f < fc; f++)
	{
		d2Vector doseVec0 = Damage::findDoseWeight(
				f, refFrame, tileSize, powSum, Nk, r0dmg_px, 
				0.0, 
				200.0, 
				10.0);
		
		d2Vector doseVec1  = Damage::findDoseWeight(
				f, refFrame, tileSize, powSum, Nk, r0dmg_px, 
				doseVec0[0] - 10.0, 
				doseVec0[0] + 10.0, 
				1.0);
		
		d2Vector doseVec2  = Damage::findDoseWeight(
				f, refFrame, tileSize, powSum, Nk, r0dmg_px, 
				doseVec1[0] - 1.0, 
				doseVec1[0] + 1.0, 
				0.1);
		
		d2Vector doseVec3  = Damage::findDoseWeight(
				f, refFrame, tileSize, powSum, Nk, r0dmg_px, 
				doseVec2[0] - 0.1, 
				doseVec2[0] + 0.1, 
				0.01);
				
		
		const double bestN = doseVec3[0];
		const double bestAf = doseVec3[1];
		
		if (f % 10 == 0) std::cout << std::endl << "      ";
		
		std::cout << bestN; 
		
		if (f == fc) std::cout << std::endl;
		else std::cout << ", ";
		
		frameDose[f] = bestN;
		
		if (diag)
		{		
			std::stringstream sts;
			sts << f;
			std::ofstream ofs("debug/dose-fit_"+sts.str()+".dat");
			
			for (int r = r0dmg_px; r < sh; r++)
			{
				const double est = bestAf * exp(-0.5*bestN/Nk[r]) * powSum[refFrame][r];			
				ofs << r << " " << est << "\n";
			}
			
			ofs << "\n";
			
			for (int r = r0dmg_px; r < sh; r++)
			{			
				ofs << r << " " << powSum[f][r] << "\n";
			}
		}
	}
	
	std::cout << std::endl;
	
	return frameDose;
}

d2Vector Damage::findDoseWeight(
	int f, int refFrame, int tileSize, 
	const std::vector<std::vector<double>>& powSum,
	std::vector<double>& Nk, double k0,
	double N0, double N1, double dN)
{
	if (f == refFrame) 
	{
		return d2Vector(0.0, 1.0);
	}
	
	double bestN(0.0), bestAf(0.0);
	double leastCost = std::numeric_limits<double>::max();
	
	const int s = tileSize;
	const int sh = s/2 + 1;
	
	std::vector<double> wgk(sh);
			
	for (double N = N0; N < N1; N += dN)
	{
		for (int r = 0; r < sh; r++)
		{
			wgk[r] = exp(-0.5*N/Nk[r]) * powSum[refFrame][r];
		}
		
		double num(0.0);
		double denom(0.0);
		
		// assume:
		// powSum[f][r] = af exp(-N/2Nk) powSum[refFrame][r]
		// =>
		// af = Sum_r{ powSum[f][r] exp(-N/2Nk) powSum[refFrame][r]}
		//    / Sum_r{ (exp(-N/2Nk) powSum[refFrame][r])^2 }
		
		for (int r = k0; r < sh; r++)
		{
			num += powSum[f][r] * wgk[r] * r;
			denom += wgk[r] * wgk[r] * r;
		}
		
		const double af = num / denom;
		
		double cost(0.0);
		
		for (int r = k0; r < sh; r++)
		{
			const double e = powSum[f][r] - af * wgk[r];
			cost += r * e*e;
		}
		
		if (cost < leastCost)
		{
			bestN = N;
			bestAf = af;
			leastCost = cost;
		}
	}
	
	return d2Vector(bestN, bestAf);
}



std::pair<std::vector<double>, std::vector<double>> 
Damage::fitBkFactors(
		const BufferedImage<double>& fcc3, 
		double boxSize, 
		double pixelSize, 
		int k0, int k1,
		BufferedImage<double>& plot,
		bool useL2)
{
    const int fc = fcc3.ydim;
	
	const double B1 = useL2? 50.0 : 500.0;	

	const int kc = fcc3.xdim;
	plot = BufferedImage<double>(kc,fc);
	plot.fill(0.0);
	
    std::vector<double> B_f(fc), alpha(fc);

	const double angbox = pixelSize * boxSize;
	
	for (int f = 0; f < fc; f++)
	{
		double B_k5 = sigma_to_B(5, angbox);
		
		int k0_curr = k0;
		
		d2Vector Bk;
		
		if (useL2)
		{
			Bk = findBkRec_L2(fcc3, f, boxSize, pixelSize, k0_curr, k1, B_k5, B1, 20, 4, 0.1);
		
			while ((Bk[1] < 0.5) || (Bk[0] < 0.0 && B_to_sigma(Bk[0], angbox) < k0_curr))
			{
				k0_curr--;
				
				if (k0_curr < 10)
				{
					Bk = d2Vector(-100.0, 0.0);
					break;
				}
						
				Bk = findBkRec_L2(fcc3, f, boxSize, pixelSize, k0_curr, k1, B_k5, B1, 20, 4, 0.1);
			}
		}
		else
		{
			Bk = findBkRec_FCC(fcc3, f, boxSize, pixelSize, k0_curr, k1, B_k5, B1, 20, 4, 0.1);
		
			while ((Bk[1] < 1e-2) || (Bk[0] < 0.0 && B_to_sigma(Bk[0], angbox) < k0_curr))
			{
				k0_curr--;
				
				if (k0_curr < 10)
				{
					Bk = d2Vector(-100.0, 0.0);
					break;
				}
						
				Bk = findBkRec_FCC(fcc3, f, boxSize, pixelSize, k0_curr, k1, B_k5, B1, 20, 4, 0.1);
			}
		}
			   
			
		B_f[f] = Bk[0];
		alpha[f] = Bk[1];
		
		if (useL2)
		{
			for (int x = k0_curr; x < kc; x++)
			{
				const double kA = x / (boxSize * pixelSize);
				plot(x,f) = alpha[f] * exp(4.0 * B_f[f] * kA * kA) * fcc3(x,f,2);
			}
		}
		else
		{
			for (int x = k0_curr; x < kc; x++)
			{
				const double kA = x / (boxSize * pixelSize);
				plot(x,f) = alpha[f] * exp(4.0 * B_f[f] * kA * kA);
			}
		}
	}

    return std::make_pair(B_f,alpha);
}

BufferedImage<double> Damage::renderFit(
		const std::vector<double>& B, 
		const std::vector<double>& scale, 
		double boxSize, double pixelSize,
		int kc)
{
	const int fc = B.size();
	
    BufferedImage<double> out(kc,fc);
	
	const double p2a = 1.0 / (boxSize * pixelSize);

    for (int k = 0; k < kc; k++)
    for (int f = 0; f < fc; f++)
    {
		const double kA = k * p2a;
		
        out(k,f) = scale[f] * exp(4.0 * B[f] * kA * kA);
    }

	return out;
}

BufferedImage<double> Damage::renderFitNrm(
		const std::vector<double> &B, 
		const std::vector<double> &scale, 
		double boxSize, double pixelSize, 
		const BufferedImage<double> &FCC3)
{
	const int kc = FCC3.xdim;
	const int fc = FCC3.ydim;
	
    BufferedImage<double> out(kc,fc);
	
	const double p2a = 1.0 / (boxSize * pixelSize);

    for (int k = 0; k < kc; k++)
    for (int f = 0; f < fc; f++)
    {
		const double kA = k * p2a;
		
        out(k,f) = scale[f] * exp(4.0 * B[f] * kA * kA) * FCC3(k,f,2);
    }

	return out;
}

d2Vector Damage::findBkRec_L2(
        const BufferedImage<double>& fcc3, int f,
		double boxSize, double pixelSize,
        int k0, int k1,
        double B0, double B1,
        int steps, int depth, double q)
{
	/*for (int i = 0; i < depth; i++)
	{
		std::cout << "  ";
	}
	
	std::cout << B0 << " ... " << B1 << std::endl;*/
			
    double minErr = std::numeric_limits<double>::max();
    double bestB = B0;
    double bestScale = 1.0;

	const double p2a = 1.0 / (boxSize * pixelSize);

    for (int s = 0; s < steps; s++)
    {
        const double B = B0 + s*(B1 - B0)/(steps-1);

        // find scale
        double num = 0.0, denom = 0.0;

        for (int k = k0; k < k1; k++)
        {
			const double k_Ang = p2a * k;
			
			double p = fcc3(k,f,0);
			double g = fcc3(k,f,1) / k;
			double q = fcc3(k,f,2) * exp(4 * B * k_Ang * k_Ang);

            num += q*p / g;
            denom += q*q / g;
        }

        const double eps = 1e-20;
        double scale = denom > eps? num / denom : num / eps;
		
		if (scale < 0.0) scale = 0.0;

        double sum = 0.0;

        for (int k = k0; k < k1; k++)
        {
			const double k_Ang = p2a * k;
			
			double p = fcc3(k,f,0);
			double g = fcc3(k,f,1) / k;
			double q = fcc3(k,f,2) * scale * exp(4 * B * k_Ang * k_Ang);
			
            const double d = q - p;
            sum += d * d / g;
        }

        if (sum < minErr)
        {
            minErr = sum;
            bestB = B;
            bestScale = scale;
        }
    }

    if (depth > 0)
    {
        const double hrange = 0.5 * (B1 - B0);
        double nextB0 = bestB - q*hrange;
        double nextB1 = bestB + q*hrange;
		
		if (nextB0 < B0) nextB0 = B0;
		if (nextB1 > B1) nextB1 = B1;

        return findBkRec_L2(fcc3, f, boxSize, pixelSize, k0, k1, nextB0, nextB1, steps, depth - 1, q);
    }

	return d2Vector(bestB, bestScale);
}

d2Vector Damage::findBkRec_FCC(
        const BufferedImage<double>& fcc3, int f,
		double boxSize, double pixelSize,
        int k0, int k1,
        double B0, double B1,
        int steps, int depth, double q)
{
	/*for (int i = 0; i < depth; i++)
	{
		std::cout << "  ";
	}
	
	std::cout << B0 << " ... " << B1 << std::endl;*/
			
    double minErr = std::numeric_limits<double>::max();
    double bestB = B0;
    double bestScale = 1.0;

	const double p2a = 1.0 / (boxSize * pixelSize);

    for (int s = 0; s < steps; s++)
    {
        const double B = B0 + s*(B1 - B0)/(steps-1);

        // find scale
        double num = 0.0, denom = 0.0;

        for (int k = k0; k < k1; k++)
        {
			const double k_Ang = p2a * k;
			
			const double wg = sqrt(fcc3(k,f,1)*fcc3(k,f,2));
			
			if (wg > 0.0)
			{
				const double p = fcc3(k,f,0) / wg;				
				const double q = exp(4 * B * k_Ang * k_Ang);
			
				num += q*p;
				denom += q*q;
			}
        }

        const double eps = 1e-20;
        double scale = denom > eps? num / denom : num / eps;
		
		if (scale < 0.0) scale = 0.0;

        double sum = 0.0;

        for (int k = k0; k < k1; k++)
        {
			const double k_Ang = p2a * k;			
			const double wg = sqrt(fcc3(k,f,1)*fcc3(k,f,2));
			
			if (wg > 0.0)
			{
				const double p = fcc3(k,f,0) / wg;				
				const double q = scale * exp(4 * B * k_Ang * k_Ang);
			
				const double d = q - p;
				sum += d * d;
			}
        }

        if (sum < minErr)
        {
            minErr = sum;
            bestB = B;
            bestScale = scale;
        }
    }

    if (depth > 0)
    {
        const double hrange = 0.5 * (B1 - B0);
        double nextB0 = bestB - q*hrange;
        double nextB1 = bestB + q*hrange;
		
		if (nextB0 < B0) nextB0 = B0;
		if (nextB1 > B1) nextB1 = B1;

        return findBkRec_FCC(fcc3, f, boxSize, pixelSize, k0, k1, nextB0, nextB1, steps, depth - 1, q);
    }

	return d2Vector(bestB, bestScale);
}

double Damage::B_to_sigma(double B, double angbox)
{
	/* B =  -angbox² / (8 * sigma²);
		<=>
	   sigma = sqrt( -angbox² / (8B) )	
	  
	*/
	
	return sqrt( - angbox * angbox / (8.0 * B) );
}

double Damage::sigma_to_B(double sigma, double angbox)
{
	/*
	  4Bx² = -k² / (2 sig²)
	  x = k / (pixelSize * boxSize)
	  => 
	  B = -k² / (8 sig² x²) = -(pixelSize * boxSize)² / (8 sig²)
	  
	  */
	
	return - angbox * angbox / (8.0 * sigma * sigma);
}

BufferedImage<float> Damage::computeBfactorWeights(
        const std::vector<double>& bFacs, 
		const std::vector<double>& aFacs, 
		int boxSize, int pixelSize, 
		bool normalize)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	
    const int fc = bFacs.size();

    BufferedImage<float> out(sh,s,fc);
	
	const double p2a = 1.0 / (pixelSize * boxSize);
	const double a = 4.0 * p2a * p2a;

    for (int f = 0; f < fc; f++)
    {
        for (int y = 0; y < s; y++)
        for (int x = 0; x < sh; x++)
        {
            double yy = y < s/2? y : y - s;
            double r2 = x*x + yy*yy;
			
			out(x,y,f) = aFacs[f] * exp(a * bFacs[f] * r2);
        }
    }

    if (normalize)
    {
        for (int y = 0; y < s; y++)
        for (int x = 0; x < sh; x++)
        {
            double sum = 0.0;

            for (int f = 0; f < fc; f++)
            {
                sum += out(x,y,f);
            }

            for (int f = 0; f < fc; f++)
            {
                out(x,y,f) /= sum;
            }
        }
    }

	return out;
}

void Damage::renormalise(
		std::vector<std::vector<double>>& B_t, 
		std::vector<std::vector<double>>& k_t)
{
	const int tc = B_t.size();
	
	double B_max = -std::numeric_limits<double>::max();
	
	for (int t = 0; t < tc; t++)
	{
		const int fc = B_t[t].size();
		
		for (int f = 0; f < fc; f++)
		{
			const double B = B_t[t][f];
			
			if (k_t[t][f] > 0.0 && B > B_max) B_max = B;
		}
	}
	
	std::cout << "B_max = " << B_max << std::endl;
	
	for (int t = 0; t < tc; t++)
	{
		const int fc = B_t[t].size();
		
		for (int f = 0; f < fc; f++)
		{
			B_t[t][f] -= B_max;
		}
	}
}
