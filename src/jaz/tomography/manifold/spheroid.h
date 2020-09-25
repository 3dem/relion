#ifndef SPHEROID_H
#define SPHEROID_H

#include "manifold.h"

class Spheroid : public Manifold
{
	public:
		
		
		Spheroid(const std::vector<double>& parameters, int index);
		
		
			gravis::d3Vector position;
			std::vector<double> SH_coefficients;
			
		
		std::vector<RigidAlignment> sampleParticles(
										double spacing,
										double depth,
										double minTilt,
										double maxTilt);
		
		std::vector<double> getParameters() const;
		
		static std::string getTypeName();
};

#endif
