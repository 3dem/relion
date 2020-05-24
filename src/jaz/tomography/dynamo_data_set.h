#ifndef DYNAMO_DATA_SET_H
#define DYNAMO_DATA_SET_H

#include <src/jaz/tomography/dynamo/dynamo_particle.h>


class DynamoDataSet
{
	public:
		
		DynamoDataSet();
		DynamoDataSet(std::string filename);
		
			std::vector<DynamoParticle> particles;
	
		
		std::vector<std::vector<int>> splitByTomogram() const;
		
		int getTotalParticleNumber() const;
		gravis::d3Vector getPosition(long int particle_id) const;
		gravis::d3Matrix getMatrix3x3(long int particle_id) const;
		gravis::d4Matrix getMatrix4x4(long int particle_id, double w, double h, double d) const;	
		std::string getName(long int particle_id) const;
		int getHalfSet(long int particle_id) const;
		
		void moveParticleTo(long int particle_id, gravis::d3Vector pos);
		void shiftParticleBy(long int particle_id, gravis::d3Vector shift);
		
		void write(std::string fn) const;
};

#endif
