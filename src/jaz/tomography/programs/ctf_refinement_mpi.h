#ifndef RELION_CTF_REFINEMENT_MPI_H
#define RELION_CTF_REFINEMENT_MPI_H

#include "ctf_refinement.h"
#include <src/mpi.h>

class CtfRefinementProgramMpi : public CtfRefinementProgram
{
private:
	MpiNode *node;

public:

	CtfRefinementProgramMpi(int argc, char *argv[]);

	/** Destructor, calls MPI_Finalize */
	~CtfRefinementProgramMpi()
	{
		delete node;
	}

	int rank, nodeCount;

	void run();
};


#endif //RELION_CTF_REFINEMENT_MPI_H
