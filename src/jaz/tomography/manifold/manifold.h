#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <src/jaz/math/rigid_alignment.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/metadata_table.h>
#include <vector>

class Manifold
{
	public:
		
		Manifold(int index, const std::string& type) 
		:	index(index), 
			type(type) 
		{}

		virtual ~Manifold() {}
		
			int index;
			std::string type;
		
		virtual std::vector<RigidAlignment> sampleParticles(
					double spacing,
					double depth,
					double minTilt,
					double maxTilt,
					bool sample_present_wedge,
					bool sample_missing_wedge,
					const std::vector<gravis::d4Matrix>& projections) const = 0;
			
		virtual std::vector<double> getParameters() const = 0;

};

#endif
