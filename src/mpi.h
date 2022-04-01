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
/***************************************************************************
 * Authors: J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef MPI_H_
#define MPI_H_
#include <mpi.h>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "src/error.h"
#include "src/macros.h"

#define MPITAG_JOB_REQUEST 0
#define MPITAG_JOB_REPLY 1
#define MPITAG_METADATA 2
#define MPITAG_RANDOMSEED 3
#define MPITAG_IMAGE 4
#define MPITAG_IMAGE_SIZE 41
#define MPITAG_PACK 5
#define MPITAG_RFLOAT 6
#define MPITAG_INT 7
#define MPITAG_IDENTIFIER 8
#define MPITAG_BCAST 9
#define MPITAG_WAIT 10
#define MPITAG_STOP 11

/** Class to wrapp some MPI common calls in an work node.
*
*/
class MpiNode
{

public:
	int rank, size;

	MPI_Group worldG, followerG; // groups of ranks (in practice only used to create communicators)
	MPI_Comm worldC, followerC; // communicators
	int followerRank; // index of follower within the follower-group (and communicator)

	MpiNode(int &argc, char ** argv);

	~MpiNode();

	// Only true if rank == 0
	bool isLeader() const;

	// Prints the random subset for this rank
	int myRandomSubset() const;

	// Returns the name of the host this rank is running on
	std::string getHostName() const;

	/** Wait on a barrier for the other MPI nodes */
	void barrierWait();

	/** Build in some better error handling and (hopefully better) robustness to communication failures in the MPI_Send/MPI_Recv pairs....
	 */
	int relion_MPI_Send(void *buf, std::ptrdiff_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

	int relion_MPI_Recv(void *buf, std::ptrdiff_t count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status &status);

	int relion_MPI_Bcast(void *buffer, long int count, MPI_Datatype datatype, int root, MPI_Comm comm);

	/* Better error handling of MPI error messages */
	void report_MPI_ERROR(int error_code);

};

// General function to print machinenames on all MPI nodes
void printMpiNodesMachineNames(MpiNode &node, int nthreads = 1);
#endif /* MPI_H_ */
