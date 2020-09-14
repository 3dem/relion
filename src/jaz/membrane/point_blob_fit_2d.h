#ifndef POINT_BLOB_FIT_H
#define POINT_BLOB_FIT_H

#include <src/jaz/optimization/optimization.h>
#include "blob_2d.h"

#include <omp.h>


class PointBlobFit2D : public Optimization
{
	public:
		
		PointBlobFit2D(
			const std::vector<gravis::d2Vector>& allPoints,
			gravis::d2Vector origin,
			double initial_radius, 
			double tolerance,
		    double tethering);
		
		
		
			gravis::d2Vector origin;
			double initial_radius, tolerance, tethering;
			std::vector<gravis::d2Vector> points;
			
		
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
