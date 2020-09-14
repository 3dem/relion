#ifndef TOPAZ_HELPER_H
#define TOPAZ_HELPER_H

#include <string>
#include <vector>
#include <map>
#include <src/jaz/gravis/t2Vector.h>

class TopazHelper
{
	public:
		
		struct Particle
		{
			gravis::i2Vector coordinates;
			double score;
		};
		
		static std::map<std::string, std::vector<Particle>> read(
		        std::string fileName, 
		        double minScore);
};

typedef std::map<std::string, std::vector<TopazHelper::Particle>> TopazParticleMap;

#endif
