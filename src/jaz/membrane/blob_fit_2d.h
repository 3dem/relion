#ifndef BLOB_FIT_3D_H
#define BLOB_FIT_3D_H

#include <src/jaz/optimization/optimization.h>
#include "blob_2d.h"

#include <omp.h>



class BlobFit2D : public Optimization
{
	public:
		
		BlobFit2D(
			const RawImage<float>& image,
			gravis::d2Vector position,
			double smoothingRadius,
			double priorSigma,
			double roundedness,
			int num_threads);


				const RawImage<float>& image;
				gravis::d2Vector initialPos;
				int smoothingRadius;
				BufferedImage<float> weight;
				double priorSigma2, roundedness;
				int num_threads;
			

		double f(const std::vector<double>& x, void* tempStorage) const;
		
		void computeWeight(
		        const Blob2D& blob,
		        double minRadius,
		        double maxRadius,
		        const RawImage<float>& mask);
		
		void* allocateTempStorage() const {return 0;}
		void deallocateTempStorage(void* ts) const {}
		
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
