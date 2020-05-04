#ifndef MOTIF_REFINEMENT_H
#define MOTIF_REFINEMENT_H

#include "motif_detection.h"
#include <src/jaz/optimization/optimization.h>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>

class Motif;
class PointCloud;

class MotifRefinement : public DifferentiableOptimization
{
	public:
		
		MotifRefinement(
				const Motif& motif,
				const std::vector<PointCloud>& clouds, 
				const std::vector<MotifDetection>& detections, 
				double tolerance,
				double reg);
		
		
			int mc, dc;
			
			std::vector<gravis::d3Vector> initialMotifPts, targetPts;
			std::vector<size_t> neighborInd, neighborNum, detectionInds, tomoInds;
			double tolerance, reg;
			
			std::vector<gravis::d3Vector> initialCenters, initialAngles;
					
			
		double f(const std::vector<double>& x, void* tempStorage) const;
		
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
		
		
		std::vector<double> getInitialParams();
		std::pair<Motif, std::vector<MotifDetection>> getDetections(const std::vector<double>& x);
				
};

#endif
