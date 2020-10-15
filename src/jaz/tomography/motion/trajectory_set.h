#ifndef TRAJECTORY_SET_H
#define TRAJECTORY_SET_H

#include <src/jaz/gravis/t3Vector.h>
#include <vector>


class TrajectorySet
{
	public:
		
		TrajectorySet();
		TrajectorySet(const std::vector<std::vector<gravis::d3Vector>>& shiftsInPix, double pixelSize);
		TrajectorySet(std::string filename);
				
			std::vector<std::vector<gravis::d3Vector>> shifts;
						
		void write(std::string filename);
};

#endif
