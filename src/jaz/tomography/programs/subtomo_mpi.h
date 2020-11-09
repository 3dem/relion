#ifndef RELION_SUBTOMO_MPI_H
#define RELION_SUBTOMO_MPI_H

#include <src/mpi.h>
#include "subtomo.h"

class SubtomoProgramMpi : public SubtomoProgram
{
private:
	MpiNode *node;

public:

	/** Destructor, calls MPI_Finalize */
	~SubtomoProgramMpi()
	{
		delete node;
	}

	int rank, nodeCount, verb;

	void readParameters(int argc, char *argv[]);
	void run();
};

#endif //RELION_SUBTOMO_MPI_H
