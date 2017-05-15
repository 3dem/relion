#ifndef LOCAL_SYMMETRY_MPI_H_
#define LOCAL_SYMMETRY_MPI_H_

#include "src/mpi.h"
#include "src/local_symmetry.h"
#include "src/parallel.h"

//#define DEBUG

#define MPITAG_LOCALSYM_SAMPLINGS_PACK 1

class local_symmetry_parameters_mpi: public local_symmetry_parameters
{
private:
	MpiNode *node;

public:
	/** Destructor, calls MPI_Finalize */
    ~local_symmetry_parameters_mpi()
    {
        delete node;
    }

    /** Read
     * This could take care of mpi-parallelisation-dependent variables
     */
    void read(int argc, char **argv);

    // Parallelized run function
    void run();

};

#endif /* LOCAL_SYMMETRY_MPI_H_ */
