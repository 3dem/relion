#ifndef RELION_ALIGN_MPI_H
#define RELION_ALIGN_MPI_H

#include <src/mpi.h>
#include "align.h"


class AlignProgramMpi : public AlignProgram
{
private:
	MpiNode *node;

public:

	AlignProgramMpi(int argc, char *argv[]);

	/** Destructor, calls MPI_Finalize */
	~AlignProgramMpi()
	{
		delete node;
	}

	int rank, nodeCount;

	void run();
};

#endif //RELION_ALIGN_MPI_H
