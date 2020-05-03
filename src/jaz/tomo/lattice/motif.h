#ifndef LATTICE_MOTIF_H
#define LATTICE_MOTIF_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t4Matrix.h>
#include "point_cloud.h"
#include "motif_detection.h"

class Motif
{
	public:
		
		
		static std::vector<std::pair<Motif, double>> discover(
				const std::vector<PointCloud>& clouds, 
				int neighbors, 
				int motifsNum = 10000, 
				int evalsNum = -1,
				double tolerance = 0.25,
                int number = 10,
                bool refine = true,
				int num_threads = 1);
		
		static std::pair<gravis::d4Matrix, double> alignPointsRansac(
				const std::vector<gravis::d3Vector>& pts_a,
				const std::vector<gravis::d3Vector>& pts_b,
                double tolerance, double offset = 1.0);

        static std::pair<gravis::d4Matrix, double> refineAlignment(
                const std::vector<gravis::d3Vector>& pts_a,
                const std::vector<gravis::d3Vector>& pts_b,
                double tolerance,
                const gravis::d4Matrix initialRot);

        static double evaluateAlignment(
                const std::vector<gravis::d3Vector>& pts_a,
                const std::vector<gravis::d3Vector>& pts_b,
                double tolerance,
                const gravis::d4Matrix A);
		
		
		
		
		Motif();		
		Motif(const std::vector<gravis::d3Vector>& points);
		
		
			std::vector<gravis::d3Vector> points;
			
			
		std::vector<MotifDetection> detectIn(
				const std::vector<PointCloud>& clouds, 
				double minScore,
                double tolerance,
                bool refine = true) const;
		
		std::vector<Mesh> visualize(
				const std::vector<MotifDetection>& occurences) const;
		
		Mesh toMesh(
				double length = 0.42, 
				double pointedness = 0.2) const;		
		
		std::pair<Motif, std::vector<MotifDetection>> refine(
				const std::vector<PointCloud>& clouds,
				const std::vector<MotifDetection>& detections,
				double tolerance,
				double reg);
				
};


#endif
