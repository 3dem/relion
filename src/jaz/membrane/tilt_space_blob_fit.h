#ifndef TILT_SPACE_BLOB_FIT_H
#define TILT_SPACE_BLOB_FIT_H

#include <src/jaz/optimization/optimization.h>
#include "blob_3d.h"

class TiltSpaceBlobFit : public DifferentiableOptimization
{
	public:
		
		TiltSpaceBlobFit(
		        int sh_bands,
		        const RawImage<float>& correlation,
		        const RawImage<gravis::d3Vector>& directions_xz );
		
		
			int sh_bands;
			const RawImage<float>& correlation;		
			BufferedImage<double> basis;
			
		
		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
		
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
