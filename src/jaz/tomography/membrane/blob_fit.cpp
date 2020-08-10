#include "blob_fit.h"

BlobFit::BlobFit(
		const Tomogram& tomogram,
		const std::vector<gravis::d4Vector>& positions,
		int id, double mean_radius,
		int outer_radius, int sh_bands, bool useMasks, double priorSigma,
		int num_threads)
:	tomogram(tomogram),
	initialPos(positions[id].x, positions[id].y, positions[id].z),
	outer_radius(outer_radius), 
	sh_bands(sh_bands),
	useMasks(useMasks),
	num_threads(num_threads),
	priorSigma2(priorSigma*priorSigma)
{
	const int fc = tomogram.stack.zdim;
	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;
		
	const double maskedVal = 1e-4;
	
	masks.resize(fc);
			
	gravis::d4Vector pw(initialPos);
	
	for (int f = 0; f < fc; f++)
	{
		if (!useMasks) continue;
		
		masks[f] = BufferedImage<float>(w, h);
		masks[f].fill(1.f);
		
		for (int i = 0; i < positions.size(); i++)
		{
			if (i == id) continue;
			
			gravis::d4Vector piw(positions[id].x, positions[id].y, positions[id].z, 1.0);
			gravis::d4Vector pii = tomogram.projectionMatrices[f] * piw;
			
			for (int y = pii.y - mean_radius; y <= pii.y + mean_radius; y++)
			for (int x = pii.x - mean_radius; x <= pii.x + mean_radius; x++)
			{
				if (x < 0 || x >= w || y < 0 || y >= h) continue;
				
				double dx = x - pii.x;
				double dy = y - pii.y;
				
				if (dx*dx + dy*dy < mean_radius*mean_radius)
				{
					masks[f](x,y) = maskedVal;
				}
			}
			
		}
	}

	std::cout << "initialPos = " << initialPos << std::endl;
}

double BlobFit::f(const std::vector<double>& x, void* tempStorage) const
{
	SphericalHarmonics* sh = (SphericalHarmonics*) tempStorage;
	Blob blob(x, outer_radius, sh);

	const int fc = tomogram.stack.zdim;
	
	int padding = 4069;
	std::vector<double> per_thread(padding * num_threads, 0.0);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		std::vector<double> radAvg = blob.radialAverage(
				tomogram, f, -1, useMasks? &masks[f] : 0);

		per_thread[padding*t] += blob.radialAverageError(
				tomogram, f, radAvg, useMasks? &masks[f] : 0);
	}
	
	double out = 0.0;
	
	for (int t = 0; t < num_threads; t++)
	{
		out += per_thread[padding*t];
	}
	
	out += (blob.center - initialPos).norm2() / priorSigma2;
	
	return out;
}

void BlobFit::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	SphericalHarmonics* sh = (SphericalHarmonics*) tempStorage;
	Blob blob(x, outer_radius, sh);

	const int fc = tomogram.stack.zdim;
	const int cc = x.size();
	
	for (int i = 0; i < cc; i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int f = 0; f < fc; f++)
	{
		std::vector<double> radAvg = blob.radialAverage(tomogram, f);
		std::vector<double> grad = blob.radialAverageErrorGrad(tomogram, f, radAvg);

		for (int i = 0; i < cc; i++)
		{
			gradDest[i] += grad[i];
		}
	}
}
