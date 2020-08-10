#ifndef BLOB_FIT_H
#define BLOB_FIT_H

#include <src/jaz/optimization/optimization.h>
#include "blob.h"

#include <omp.h>

#include <src/spherical-harmonics/SphericalHarmonics.h>

using namespace au::edu::anu::qm::ro;


class BlobFit : public DifferentiableOptimization
{
	public:
		
		BlobFit(const Tomogram& tomogram,
				const std::vector<gravis::d4Vector>& positions,
				int id, double mean_radius,
				int outer_radius, int sh_bands, bool useMasks, double priorSigma,
				int num_threads);
		
			const Tomogram& tomogram;
			gravis::d3Vector initialPos;
			int outer_radius, sh_bands;
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

#endif
