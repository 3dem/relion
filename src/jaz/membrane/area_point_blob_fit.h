#ifndef AREA_POINT_BLOB_FIT_H
#define AREA_POINT_BLOB_FIT_H

#include <src/jaz/optimization/optimization.h>
#include "blob_2d.h"

#include <omp.h>


class AreaPointBlobFit : public Optimization
{
	public:
		
		AreaPointBlobFit(
			const RawImage<float>& pointDensity,
			const RawImage<float>& regionalTerm,
			gravis::d2Vector origin,
			double initial_radius, 
			double tolerance,
		    double tethering,
			double aspectCost,
		    double contrastCost);
		
		
		
			gravis::d2Vector origin;
			int radialSamples;
			double initial_radius, tolerance, tethering, aspectCost, contrastCost;
			const RawImage<float>& pointDensity;
			BufferedImage<float> regionalTermSquare;
			
			gravis::i2Vector boxOrigin, boxSize;
			double boxArea;
			
		
		double f(const std::vector<double>& x, void* tempStorage) const;
		
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
		
		BufferedImage<float> visualise(const std::vector<double>& x, int w, int h);
};

#endif
