#include "catalogue.h"
#include <fstream>

Catalogue::Catalogue()
{
}

Catalogue::Catalogue(std::string filename)
{
	std::ifstream ifs(filename);
	
	const int maxlen = 4096;
	char buffer[maxlen];
	
	while (ifs.getline(buffer, maxlen))
	{
		particles.push_back(DynamoParticle(buffer));
	}
}

std::vector<std::vector<DynamoParticle>> Catalogue::splitByTomogram()
{
	int maxTomo = 0;
	
	for (int p = 0; p < particles.size(); p++)
	{
		const int t = particles[p].tomo;
		
		if (t > maxTomo) maxTomo = t;
	}
	
	std::vector<std::vector<DynamoParticle>> out(maxTomo+1);
	
	for (int p = 0; p < particles.size(); p++)
	{
		out[particles[p].tomo].push_back(particles[p]);
	}
	
    return out;
}

void Catalogue::write(std::string fn)
{
    std::ofstream ofs(fn);

    for (int p = 0; p < particles.size(); p++)
    {
        ofs << particles[p] << '\n';
    }
}
