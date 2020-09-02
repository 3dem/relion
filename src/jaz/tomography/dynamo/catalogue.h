#ifndef CATALOGUE_H
#define CATALOGUE_H

#include "dynamo_particle.h"

class Catalogue
{
	public:

        Catalogue();
        Catalogue(std::string filename);

			std::vector<DynamoParticle> particles;
		
		std::vector<std::vector<DynamoParticle>> splitByTomogram();
        void write(std::string fn);
};

#endif
