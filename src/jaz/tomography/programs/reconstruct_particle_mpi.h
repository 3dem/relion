#ifndef RECONSTRUCT_PARTICLE_MPI_H
#define RECONSTRUCT_PARTICLE_MPI_H

#include <src/mpi.h>
#include "reconstruct_particle.h"

class ReconstructParticleProgramMpi : public ReconstructParticleProgram
{
	private:
	MpiNode *node;

	public:

		//ReconstructParticleProgramMpi(int rank, int nodeCount);

		/** Destructor, calls MPI_Finalize */
		~ReconstructParticleProgramMpi()
		{
			delete node;
		}

		int rank, nodeCount, verb;
		std::string tmpOutRootBase;

		void readParameters(int argc, char *argv[]);
		void run();
};

#endif
