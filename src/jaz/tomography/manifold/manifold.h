#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <src/jaz/math/rigid_alignment.h>
#include <src/metadata_table.h>
#include <vector>

class Manifold
{
	public:
		
		Manifold(int index, const std::string& type) 
		:	index(index), 
			type(type) 
		{}
		
			int index;
			std::string type;
		
		virtual std::vector<RigidAlignment> sampleParticles(
					double spacing,
					double depth,
					double minTilt,
					double maxTilt) const = 0;
			
		virtual std::vector<double> getParameters() const = 0;

};

#endif
