/*
 * autopicker_mpi.h
 *
 *  Created on: Sep 18, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#ifndef AMYLOID_FINDER_MPI_H_
#define AMYLOID_FINDER_MPI_H_

#include "src/mpi.h"
#include "src/amyloid_finder.h"
#include "src/parallel.h"

class AmyloidFinderMpi: public AmyloidFinder
{
private:
	MpiNode *node;

public:
	/** Destructor, calls MPI_Finalize */
	~AmyloidFinderMpi()
	{
		delete node;
	}

	/** Read
	 * This could take care of mpi-parallelisation-dependent variables
	 */
	void read(int argc, char **argv);

	// Set device-affinity
	void deviceInitialise();

	// Parallelized run function
	void run();

	int getRank()
	{
		return(node->rank);
	}

	MpiNode * getNode()
	{
		return(node);
	}
};

#endif /* AMYLOID_FINDER_MPI_H_ */
