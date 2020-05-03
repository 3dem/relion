#ifndef ASTIGMATISM_FIT_H
#define ASTIGMATISM_FIT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>


class AstigmatismFit : public Optimization
{
    public:
		
		AstigmatismFit(const std::vector<BufferedImage<float>>& spectra,
					   int num_threads, int k0, int k1,
					   double Cs_px);
		
			const std::vector<BufferedImage<float>>& spectra;
			int num_threads, k0, k1;
			double Cs_px;
		
		
		double f(const std::vector<double>& x, void* tempStorage) const;
		void report(int iteration, double cost, const std::vector<double>& x) const;
		
		
		void averageAlongIso(
				int f, 
				double alpha, double beta, 
				std::vector<double>& accum, 
				std::vector<double>& weight,
				int boundary = 2) const;
		
		double compareWithExpansion(
				int f, 
				double alpha, double beta, 
				const std::vector<double>& accum,
				int boundary = 2) const;
		
		BufferedImage<double> computeExpansion(
				int f, 
				double alpha, double beta, 
				const std::vector<double>& accum,
				int boundary = 2) const;
		
};

#endif
