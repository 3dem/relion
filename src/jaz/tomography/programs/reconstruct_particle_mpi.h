#ifndef RECONSTRUCT_PARTICLE_MPI_H
#define RECONSTRUCT_PARTICLE_MPI_H

#include "reconstruct_particle.h"

class ReconstructParticleProgramMpi : public ReconstructParticleProgram
{
	public:

		ReconstructParticleProgramMpi(int rank, int nodeCount);

		int rank, nodeCount;

		void readParameters(int argc, char *argv[]);
		void run();
};

#endif
