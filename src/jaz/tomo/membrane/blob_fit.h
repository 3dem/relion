#ifndef BLOB_FIT_H
#define BLOB_FIT_H

#include <src/jaz/optimization/optimization.h>
#include "blob.h"

#include <omp.h>

#include <src/spherical-harmonics/SphericalHarmonics.h>

using namespace au::edu::anu::qm::ro;


template <typename T>
class BlobFit : public DifferentiableOptimization
{
	public:
		
		BlobFit(const TomoStack<T>& ts, 
				const std::vector<gravis::d3Vector>& positions, 
				int id, double mean_radius,
				int outer_radius, int sh_bands, bool useMasks, double priorSigma,
				int num_threads);
		
			const TomoStack<T>& ts;
			gravis::d3Vector initialPos;
			int outer_radius, sh_bands;
			std::vector<bool> goodViews;
			std::vector<BufferedImage<float>> masks;
			bool useMasks;
			int num_threads;
			double priorSigma2;
			
		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
		
		void* allocateTempStorage() const {return new SphericalHarmonics(sh_bands);}
		void deallocateTempStorage(void* ts) const {delete (SphericalHarmonics*) ts;}
		
		void report(int iteration, double cost, const std::vector<double>& x) const 
		{
			std::cout << "it " << iteration << ", " << cost << ":\n";
			
			for (int i = 0; i < x.size(); i++)
			{
				std::cout << x[i] << "  ";
			}
			
			std::cout << std::endl;
		}
};

template <typename T>
BlobFit<T>::BlobFit(
		const TomoStack<T>& ts, 
		const std::vector<gravis::d3Vector>& positions, 
		int id, double mean_radius,
		int outer_radius, int sh_bands, bool useMasks, double priorSigma,
		int num_threads)
:	ts(ts), 
	initialPos(positions[id]), 
	outer_radius(outer_radius), 
	sh_bands(sh_bands),
	goodViews(ts.images.size()),
	useMasks(useMasks),
	num_threads(num_threads),
	priorSigma2(priorSigma*priorSigma)
{
	const int fc = ts.images.size();
	const int w = ts.images[0].xdim;
	const int h = ts.images[0].ydim;
		
	const double maskedVal = 1e-4;
	
	masks.resize(fc);
			
	gravis::d4Vector pw(initialPos);
	
	for (int f = 0; f < fc; f++)
	{
		gravis::d4Vector pi = ts.worldToImage[f] * pw;
		
		goodViews[f] = pi.x > outer_radius && pi.x + outer_radius < ts.images[f].xdim
				&& pi.y > outer_radius && pi.y + outer_radius < ts.images[f].ydim;
		
		if (!useMasks) continue;
		
		masks[f] = BufferedImage<float>(w, h);
		masks[f].fill(1.f);
		
		for (int i = 0; i < positions.size(); i++)
		{
			if (i == id) continue;
			
			gravis::d4Vector piw(positions[i]);
			gravis::d4Vector pii = ts.worldToImage[f] * piw;
			
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
}
	
template <typename T>
double BlobFit<T>::f(const std::vector<double>& x, void* tempStorage) const
{
	SphericalHarmonics* sh = (SphericalHarmonics*) tempStorage;
	Blob blob(x, outer_radius, sh);
	
	const int fc = ts.images.size();	
	
	int padding = 4069;
	std::vector<double> per_thread(padding * num_threads, 0.0);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		//std::cout << t << "/" << num_threads << std::endl;
		
		if (goodViews[f])
		{
			std::vector<double> radAvg = blob.radialAverage(ts, f, -1, useMasks? &masks[f] : 0);
			per_thread[padding*t] += blob.radialAverageError(ts, f, radAvg, useMasks? &masks[f] : 0);
		}
	}
	
	double out = 0.0;
	
	for (int t = 0; t < num_threads; t++)
	{
		out += per_thread[padding*t];
	}
	
	out += (blob.center - initialPos).norm2() / priorSigma2;
	
	return out;
}

template<typename T>
void BlobFit<T>::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	SphericalHarmonics* sh = (SphericalHarmonics*) tempStorage;
	Blob blob(x, outer_radius, sh);
	
	const int fc = ts.images.size();
	const int cc = x.size();
	
	for (int i = 0; i < cc; i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int f = 0; f < fc; f++)
	{
		if (goodViews[f])
		{
			std::vector<double> radAvg = blob.radialAverage(ts, f);
			std::vector<double> grad = blob.radialAverageErrorGrad(ts, f, radAvg);
			
			for (int i = 0; i < cc; i++)
			{
				gradDest[i] += grad[i];
			}
		}
	}
}

#endif
