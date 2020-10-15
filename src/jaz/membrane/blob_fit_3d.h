#ifndef BLOB_FIT_3D_H
#define BLOB_FIT_3D_H

#include <src/jaz/optimization/optimization.h>
#include "blob_3d.h"

#include <omp.h>

#include <src/spherical-harmonics/SphericalHarmonics.h>

using namespace au::edu::anu::qm::ro;


class BlobFit3D : public Optimization
{
	public:
		
		BlobFit3D(
			const Tomogram& tomogram,
			gravis::d3Vector position,
			int sh_bands,
			double radius,
			double thickness,
			const std::vector<gravis::d4Vector>& allSpheres,
			const std::vector<gravis::d3Vector>& fiducials,
			double fiducialRadius,
			double priorSigma,
			int num_threads);


				const Tomogram& tomogram;
				gravis::d3Vector initialPos;
				int outer_radius, sh_bands;
				BufferedImage<float> weight;
				double priorSigma2;
				int num_threads;
			

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
