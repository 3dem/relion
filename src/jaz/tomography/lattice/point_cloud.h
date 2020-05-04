#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/mesh/mesh.h>
#include <vector>

class PointCloud
{
	public:
		
		PointCloud();
		
			std::vector<gravis::d3Vector> vertices;
			std::vector<std::vector<int>> neighbors;
		
		void readCsv(std::string fn);
		void findNeighbors(double minDist, double maxDist, int expected = 6);
		std::vector<gravis::d3Vector> getNeighborhood(int vertex) const;
		Mesh visualizeNeighbors();
		
		std::vector<size_t> distanceHist(int bins, double maxDist);
};


#endif
