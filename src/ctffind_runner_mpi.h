/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef CTFFIND_RUNNER_MPI_H_
#define CTFFIND_RUNNER_MPI_H_

#include "src/mpi.h"
#include "src/ctffind_runner.h"
#include "src/parallel.h"

class CtffindRunnerMpi: public CtffindRunner
{
public:
	MpiNode *node;

	/** Destructor, calls MPI_Finalize */
    ~CtffindRunnerMpi()
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

#endif /* CTFFIND_RUNNER_MPI_H_ */
