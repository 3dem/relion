/***************************************************************************
 *
 * Author: "Jasenko Zivanov & Sjors H.W. Scheres"
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

#ifndef MOTION_REFINER_MPI_H_
#define MOTION_REFINER_MPI_H_

#include <src/mpi.h>
#include <src/parallel.h>

#include "motion_refiner.h"

class MotionRefinerMpi : public MotionRefiner
{
private:
	MpiNode *node;

public:
	/** Destructor, calls MPI_Finalize */
    ~MotionRefinerMpi()
    {
        delete node;
    }

    /** Read
     * This could take care of mpi-parallelisation-dependent variables
     */
    void read(int argc, char **argv);

	// For the first fraction, align movies and update the FCC
	void runWithFccUpdate();
	
	// For the second, estimate the B factors from the FCC and 
	// then align and immediately recombine the following micrographs
	void runWithRecombination();
	
};


#endif
