#include "topaz_helper.h"
#include <src/macros.h>

using namespace gravis;

std::map<std::string, std::vector<TopazHelper::Particle>> TopazHelper::read(
        std::string fileName, double minScore)
{
	std::ifstream file(fileName);
	
	if (!file)
	{
		REPORT_ERROR("Unable to read: " + fileName);
	}
	
	std::map< std::string, std::vector<Particle> > out;
	
	std::string line;
	
	std::string lastName = "";
	std::vector<Particle> currentParticles;
	
	// skip first line
	std::getline(file, line);
	
	while (std::getline(file, line))
	{
		std::stringstream sts;
		sts << line;
		
		std::string name;
		sts >> name;
		
		if (name != lastName)
		{
			if (currentParticles.size() > 0)
			{
				out[lastName] = currentParticles;
				currentParticles.clear();
			}
			
			lastName = name;
		}
		
		Particle p;
		
		sts >> p.coordinates.x;
		sts >> p.coordinates.y;
		sts >> p.score;
		
		if (p.score >= minScore)
		{
			currentParticles.push_back(p);
		}
	}
	
	if (currentParticles.size() > 0)
	{
		out[lastName] = currentParticles;
	}
	
	return out;
}
