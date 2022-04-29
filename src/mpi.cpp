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
 * 02111-1307 USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/mpi.h"

// maximum amount of data that can be sent in MPI
// 512 MB is already on a safe side. If this still causes OpenMPI crash,
// please try "mpirun --mca pml ob1".
const int RELION_MPI_MAX_SIZE = 512 * 1024 * 1024;

//#define MPI_DEBUG

//------------ MPI ---------------------------
MpiNode::MpiNode(int &argc, char ** argv)
{
	//MPI Initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	// Handle errors
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	// Set up Follower communicator -----------------------------------------
	MPI_Comm_group(MPI_COMM_WORLD, &worldG);
	int mstr[1] = {0};
	MPI_Group_excl(worldG, 1, mstr, &followerG); // exclude leader
	MPI_Comm_create(MPI_COMM_WORLD, followerG, &followerC);
	if (rank != 0)
	{
		MPI_Group_rank(followerG, &followerRank);
	}
	else
	{
		followerRank = -1;
	}
}

MpiNode::~MpiNode()
{
	MPI_Finalize();
}

bool MpiNode::isLeader() const
{
	return rank == 0;
}

int MpiNode::myRandomSubset() const
{
	if (rank == 0)
		return 0;
	else
		return (rank % 2 == 0) ? 2 : 1;
}

std::string MpiNode::getHostName() const
{
	char nodename[64] = "undefined";
	gethostname(nodename,sizeof(nodename));
	std::string result(nodename);
	return result;

}

void MpiNode::barrierWait()
{
  MPI_Barrier(MPI_COMM_WORLD);
}

// MPI_TEST will be executed every this many seconds: so this determines the minimum time taken for every send operation!!
//#define VERBOSE_MPISENDRECV
int MpiNode::relion_MPI_Send(void *buf, std::ptrdiff_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
	int result(0);
	RFLOAT start_time = MPI_Wtime();

//#define ONLY_NORMAL_SEND
//#ifdef ONLY_NORMAL_SEND
	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);
	const std::ptrdiff_t blocksize(RELION_MPI_MAX_SIZE);
	const std::ptrdiff_t totalsize(count * unitsize);
	if (totalsize <= blocksize )
	{
#ifdef MPI_DEBUG
		std::cout << "relion_MPI_Send: rank = " << rank << " count = " << count << " dest = " << dest << " tag = " << tag << " comm = " << comm << std::endl;
#endif
		result = MPI_Send(buf, count, datatype, dest, tag, comm);
		if (result != MPI_SUCCESS)
			report_MPI_ERROR(result);
	}
	else
	{
		char* const buffer(reinterpret_cast<char*>(buf));
		const std::ptrdiff_t ntimes(totalsize/blocksize);
		const std::ptrdiff_t nremain(totalsize%blocksize);
		std::ptrdiff_t i(0);
		for(; i < ntimes; ++i)
		{
#ifdef MPI_DEBUG
			std::cout << "relion_MPI_Send: rank = " << rank << " blocksize = " << blocksize << " dest = " << dest << " tag = " << tag << " comm = " << comm << std::endl;
#endif
			result = MPI_Send(buffer + i * blocksize, blocksize, MPI_CHAR, dest, tag, comm);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
		if(nremain > 0)
		{
#ifdef MPI_DEBUG
			std::cout << "relion_MPI_Send: rank = " << rank << " nremain = " << nremain << " dest = " << dest << " tag = " << tag << " comm = " << comm << std::endl;
#endif
			result = MPI_Send(buffer + i * blocksize, nremain, MPI_CHAR, dest, tag, comm);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
	}

#ifdef VERBOSE_MPISENDRECV
	if (count > 100)
		std::cerr <<" relion_MPI_Send: message to " << dest << " of size "<< count << " arrived in " << MPI_Wtime() - start_time << " seconds" << std::endl;
#endif
	return result;
}

int MpiNode::relion_MPI_Recv(void *buf, std::ptrdiff_t count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status &status) {
	int result;
	MPI_Request request;
	RFLOAT current_time = MPI_Wtime();
	RFLOAT start_time = current_time;

	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);
	const std::ptrdiff_t blocksize(RELION_MPI_MAX_SIZE);
	const std::ptrdiff_t totalsize(count * unitsize);
	if (totalsize <= blocksize)
	{
#ifdef MPI_DEBUG
		std::cout << "relion_MPI_Recv: rank = " << rank << " count = " << count << " source = " << source << " tag = " << tag << " comm = " << comm << std::endl;
#endif
		int result_irecv = MPI_Irecv(buf, count, datatype, source, tag, comm, &request);
		if (result_irecv != MPI_SUCCESS)
			report_MPI_ERROR(result_irecv);

		result = MPI_Wait(&request, &status);
		if (result != MPI_SUCCESS)
			report_MPI_ERROR(result);
	}
	else
	{
		char* const buffer(reinterpret_cast<char*>(buf));
		const std::ptrdiff_t ntimes(totalsize / blocksize);
		const std::ptrdiff_t nremain(totalsize % blocksize);
		std::ptrdiff_t i(0);
		for(; i < ntimes; ++i)
		{
#ifdef MPI_DEBUG
			std::cout << "relion_MPI_Recv: rank = " << rank << " blocksize = " << blocksize << " source = " << source << " tag = " << tag << " comm = " << comm << std::endl;
#endif
			int result_irecv = MPI_Irecv(buffer + i * blocksize, blocksize, MPI_CHAR, source, tag, comm, &request);
			if (result_irecv != MPI_SUCCESS)
				report_MPI_ERROR(result_irecv);

			result = MPI_Wait(&request, &status);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
		if (nremain > 0)
		{
#ifdef MPI_DEBUG
			std::cout << "relion_MPI_Recv: rank = " << rank << " nremain = " << nremain << " source = " << source << " tag = " << tag << " comm = " << comm << std::endl;
#endif
			int result_irecv = MPI_Irecv(buffer + i * blocksize, nremain, MPI_CHAR, source, tag, comm, &request);
			if (result_irecv != MPI_SUCCESS)
				report_MPI_ERROR(result_irecv);

			result = MPI_Wait(&request, &status);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
	}

#ifdef VERBOSE_MPISENDRECV
	if (count > 100)
		std::cerr <<" relion_MPI_Recv: message from "<<source << " of size "<< count <<" arrived in " << MPI_Wtime() - start_time << " seconds" << std::endl;
#endif
		return result;
}


int MpiNode::relion_MPI_Bcast(void *buffer, long int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
	int result;
	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);

	const long blocksize(RELION_MPI_MAX_SIZE);
	const long totalsize(count * unitsize);

#ifdef MPI_DEBUG
	std::cout << "relion_MPI_Bcast: rank = " << rank << " count = " << count << " root = " << root << " comm = " << comm << std::endl;
#endif
	if (count < 0)
		report_MPI_ERROR(MPI_ERR_COUNT);  // overflow
	if (totalsize <= blocksize)
	{
		result = MPI_Bcast(buffer, static_cast<int>(count), datatype, root, comm);
		if (result != MPI_SUCCESS)
			report_MPI_ERROR(result);
	}
	else
	{
		int rank_in_group = rank, size_of_group = size;
		if (comm != MPI_COMM_WORLD)
		{
			MPI_Group group_of_comm;
			MPI_Comm_group(comm, &group_of_comm);
			MPI_Group_rank(group_of_comm, &rank_in_group);
			MPI_Group_size(group_of_comm, &size_of_group);
		}

#ifdef MPI_DEBUG
		std::cout << "relion_MPI_Bcast: global_rank = " << rank << " rank_in_comm = " << rank_in_group << " size_of_group = " << size_of_group << std::endl;
#endif
		if (rank_in_group == root)
		{
			for (int dest = 0; dest < size_of_group; dest++)
			{
				if (dest != root)
				{
					result = relion_MPI_Send(buffer, count, datatype, dest, MPITAG_BCAST, comm);
					if (result != MPI_SUCCESS)
						report_MPI_ERROR(result);
				}
			}
		}
		else
		{
			MPI_Status status;
			result = relion_MPI_Recv(buffer, count, datatype, root, MPITAG_BCAST, comm, status);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
	}

	return result;
}

void MpiNode::report_MPI_ERROR(int error_code)
{
	char error_string[200];
	int length_of_error_string, error_class;
	MPI_Error_class(error_code, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	fprintf(stderr, "%3d: %s\n", rank, error_string);
	MPI_Error_string(error_code, error_string, &length_of_error_string);
	fprintf(stderr, "%3d: %s\n", rank, error_string);

	std::cerr.flush();
	REPORT_ERROR("Encountered an MPI-related error, see above. Now exiting...");
}

void printMpiNodesMachineNames(MpiNode &node, int nthreads)
{
	if (node.isLeader())
	{
		std::cout << " === RELION MPI setup ===" << std::endl;
		std::cout << " + Number of MPI processes             = " << node.size << std::endl;
		if (nthreads > 1)
		{
			std::cout << " + Number of threads per MPI process   = " << nthreads << std::endl;
		std::cout << " + Total number of threads therefore   = " << nthreads * node.size << std::endl;
		}
		std::cout << " + Leader  (0) runs on host            = " << node.getHostName() << std::endl;
		std::cout.flush();
	}
	node.barrierWait();

	for (int follower = 1; follower < node.size; follower++)
	{
		if (follower == node.rank)
		{
			std::cout << " + Follower ";
			std::cout.width(5);
			std::cout << follower;
			std::cout << " runs on host            = " << node.getHostName() << std::endl;
			std::cout.flush();
		}
		node.barrierWait();
	}

	if (node.isLeader())
	{
		std::cout << " =================" << std::endl;
	}
	std::cout.flush();

	// Try to flush all std::cout of all MPI processes before proceeding...
	sleep(1);
	node.barrierWait();
}
