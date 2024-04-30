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

#ifdef USE_MPI_COLLECTIVE
 #include <vector>
#endif
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
	MPI_Comm_create_group(MPI_COMM_WORLD, followerG, 0, &followerC);
	if (rank != 0)
	{
		MPI_Group_rank(followerG, &followerRank);
#ifdef USE_MPI_COLLECTIVE
		// Create seperate communicator and rank for split_random_halves case run
		const int myColor = followerRank % 2;
		MPI_Comm_split(followerC, myColor, followerRank, &splitC);
		MPI_Comm_rank(splitC, &splitRank);
#endif
	}
	else
	{
		// Leader does not belong to follower
		followerC = MPI_COMM_NULL;
		followerRank = -1;
#ifdef USE_MPI_COLLECTIVE
		// Leader does not belong to split_random_halves case run
		splitC = MPI_COMM_NULL;
		splitRank = -1;
#endif
	}

#ifdef USE_MPI_COLLECTIVE
	std::vector<int> odd, even;
	for (int i = 1; i < size; i++)
	{
		if (i % 2 == 0)
			even.push_back(i);
		else
			odd.push_back(i);
	}
	root_evenC = MPI_COMM_NULL;
	root_oddC = MPI_COMM_NULL;
	MPI_Group_excl(worldG, odd.size(), odd.data(), &root_evenG);
	MPI_Group_excl(worldG, even.size(), even.data(), &root_oddG);
	MPI_Comm_create_group(MPI_COMM_WORLD, root_evenG, 1, &root_evenC);
	MPI_Comm_create_group(MPI_COMM_WORLD, root_oddG, 2, &root_oddC);

	constexpr std::ptrdiff_t defPtBlock = 4ULL*1024*1024*1024;	// 4 GiB, default for MPI communication
	constexpr std::ptrdiff_t defColBlock = 64*1024*1024;		// 64 MiB, default for MPI communication
	constexpr std::ptrdiff_t maxPtBlock = 8ULL*1024*1024*1024;		// 8 GiB, max for MPI pt-2-pt. Normal limit is sizeof(variable) * INT_MAX bytes
	constexpr std::ptrdiff_t maxCollBlock = 1ULL*1024*1024*1024;	// 1 GiB, max for MPI collective. Normal limit is sizeof(variable) * INT_MAX bytes

	// User input to change default MPI communication blocking size
	char* pEnvBlock = std::getenv("MAX_MPI_BLOCK");     // For both pt-2-pt and collective
	char* pEnvP2P   = std::getenv("MAX_MPI_P2P_BLOCK"); // For pt-2-pt blocking size
	char* pEnvColl  = std::getenv("MAX_MPI_COLL_BLOCK");// For collective blocking size

	const std::ptrdiff_t lBlock = (pEnvBlock == NULL) ? 0 : std::strtoul(pEnvBlock, nullptr, 10);
	const std::ptrdiff_t lP2P   = (pEnvP2P == NULL)   ? 0 : std::strtoul(pEnvP2P, nullptr, 10);
	const std::ptrdiff_t lColl  = (pEnvColl == NULL)  ? 0 : std::strtoul(pEnvColl, nullptr, 10);

	if (lBlock > 0) // "MAX_MPI_BLOCK" has precedence if it is set
	{
		if (lBlock <= maxPtBlock)
		{
			p2p_blocksize = lBlock;
			coll_blocksize = lBlock;
		}
		else
		{
			p2p_blocksize = maxPtBlock;
			coll_blocksize = maxCollBlock;
		}
	}
	else    // "MAX_MPI_BLOCK" is not set
	{
		if (lP2P > 0)   // "MAX_MPI_P2P_BLOCK" is set
		{
			if (lP2P <= maxPtBlock)
				p2p_blocksize = lP2P;
			else
				p2p_blocksize = maxPtBlock;
		}
		else
			p2p_blocksize = defPtBlock;

		if (lColl > 0)  // "MAX_MPI_COLL_BLOCK" is set
		{
			if (lColl <= maxCollBlock)
				coll_blocksize = lColl;
			else
				coll_blocksize = maxCollBlock;
		}
		else
			coll_blocksize = defColBlock;
	}
	// Make it multiple of 8 bytes
	p2p_blocksize = (p2p_blocksize / 8ULL) * 8ULL;
	coll_blocksize = (coll_blocksize / 8ULL) * 8ULL;

	if (rank == 0)
	{
		std::cout << "Using Point-to-Point MPI Blocksize = " << p2p_blocksize << " bytes" << std::endl;
		std::cout << "Using Collective MPI Blocksize     = " << coll_blocksize << " bytes" << std::endl;
	}
#endif
}

MpiNode::~MpiNode()
{
#ifdef USE_MPI_COLLECTIVE
	MPI_Comm_free(&root_evenC);
	MPI_Comm_free(&root_oddC);
	MPI_Comm_free(&splitC);
	MPI_Group_free(&root_evenG);
	MPI_Group_free(&root_oddG);
#endif
	MPI_Comm_free(&followerC);
	MPI_Group_free(&followerG);
	MPI_Group_free(&worldG);

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

#ifdef USE_MPI_COLLECTIVE
// Prints splitRank(split_random_halves case run) for the given nrank
int MpiNode::getSplitRank(const int nrank) const
{
	if (nrank == 0)
		return -1;
	else
		return (nrank % 2 == 0) ? nrank/2-1 : nrank/2;
}
#endif

std::string MpiNode::getHostName() const
{
	char nodename[64] = "undefined";
	gethostname(nodename,sizeof(nodename));
	std::string result(nodename);
	return result;

}

void MpiNode::barrierWait(MPI_Comm comm)
{
	MPI_Barrier(comm);
}

// MPI_TEST will be executed every this many seconds: so this determines the minimum time taken for every send operation!!
//#define VERBOSE_MPISENDRECV
int MpiNode::relion_MPI_Send(void *buf, std::ptrdiff_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
	int result(0);
#ifdef VERBOSE_MPISENDRECV
	RFLOAT start_time = MPI_Wtime();
#endif

//#define ONLY_NORMAL_SEND
//#ifdef ONLY_NORMAL_SEND
	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);
#ifdef USE_MPI_COLLECTIVE
	const std::ptrdiff_t blocksize(p2p_blocksize);
#else
	const std::ptrdiff_t blocksize(RELION_MPI_MAX_SIZE);
#endif
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
#ifdef VERBOSE_MPISENDRECV
	RFLOAT current_time = MPI_Wtime();
	RFLOAT start_time = current_time;
#endif

	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);
#ifdef USE_MPI_COLLECTIVE
	const std::ptrdiff_t blocksize(p2p_blocksize);
#else
	const std::ptrdiff_t blocksize(RELION_MPI_MAX_SIZE);
#endif
	const std::ptrdiff_t totalsize(count * unitsize);
	if (totalsize <= blocksize)
	{
#ifdef MPI_DEBUG
		std::cout << "relion_MPI_Recv: rank = " << rank << " count = " << count << " source = " << source << " tag = " << tag << " comm = " << comm << std::endl;
#endif
		result = MPI_Recv(buf, count, datatype, source, tag, comm, &status);
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
			result = MPI_Recv(buffer + i * blocksize, blocksize, MPI_CHAR, source, tag, comm, &status);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
		if (nremain > 0)
		{
#ifdef MPI_DEBUG
			std::cout << "relion_MPI_Recv: rank = " << rank << " nremain = " << nremain << " source = " << source << " tag = " << tag << " comm = " << comm << std::endl;
#endif
			result = MPI_Recv(buffer + i * blocksize, nremain, MPI_CHAR, source, tag, comm, &status);
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


int MpiNode::relion_MPI_Bcast(void *buffer, std::ptrdiff_t count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
	int result;
	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);

#ifdef USE_MPI_COLLECTIVE
	const long blocksize(coll_blocksize);
#else
	const long blocksize(RELION_MPI_MAX_SIZE);
#endif
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
#ifdef USE_MPI_COLLECTIVE
		const std::ptrdiff_t blockCount = coll_blocksize / static_cast<std::ptrdiff_t>(unitsize);
		const std::ptrdiff_t ntimes(count / blockCount);
		const int nremain(count % blockCount);
		std::ptrdiff_t i(0);
		for(; i < ntimes; ++i)
		{
			const std::ptrdiff_t offset(i * coll_blocksize);
			char* const Bbuffer(reinterpret_cast<char*>(buffer) + offset);
			result = MPI_Bcast(Bbuffer, static_cast<int>(blockCount), datatype, root, comm);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
		if (nremain > 0)
		{
			const std::ptrdiff_t offset(i * coll_blocksize);
			char* const Bbuffer(reinterpret_cast<char*>(buffer) + offset);
			result = MPI_Bcast(Bbuffer, nremain, datatype, root, comm);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
 #ifdef MPI_DEBUG
		std::cout << "relion_MPI_Bcast: rank = " << rank << " count = " << count << " root = " << root << " comm = " << comm << std::endl;
 #endif
#else	// End of USE_MPI_COLLECTIVE
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
#endif	// End of the original pt-2-pt implementation of Bcast
	}

	return result;
}

#ifdef USE_MPI_COLLECTIVE
int MpiNode::relion_MPI_Allreduce(void *sendB, void *recvB, std::ptrdiff_t count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
	int result;
	int unitsize(0);
	MPI_Type_size(datatype, &unitsize);

	const std::ptrdiff_t totalsize(count * static_cast<std::ptrdiff_t>(unitsize));

	if (count < 0)
		report_MPI_ERROR(MPI_ERR_COUNT);  // overflow

	if (totalsize <= coll_blocksize)
	{
		result = MPI_Allreduce(sendB, recvB, static_cast<int>(count), datatype, op, comm);
		if (result != MPI_SUCCESS)
			report_MPI_ERROR(result);
	}
	else
	{
		const std::ptrdiff_t blockCount = coll_blocksize / static_cast<std::ptrdiff_t>(unitsize);
		const std::ptrdiff_t ntimes(count / blockCount);
		const int nremain(count % blockCount);
		std::ptrdiff_t i(0);
		for(; i < ntimes; ++i)
		{
			const std::ptrdiff_t offset(i * coll_blocksize);
			char* const Sbuf(reinterpret_cast<char*>(sendB) + offset);
			char* const Rbuf(reinterpret_cast<char*>(recvB) + offset);
			result = MPI_Allreduce(Sbuf, Rbuf, static_cast<int>(blockCount), datatype, op, comm);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
		if (nremain > 0)
		{
			const std::ptrdiff_t offset(i * coll_blocksize);
			char* const Sbuf(reinterpret_cast<char*>(sendB) + offset);
			char* const Rbuf(reinterpret_cast<char*>(recvB) + offset);
			result = MPI_Allreduce(Sbuf, Rbuf, nremain, datatype, op, comm);
			if (result != MPI_SUCCESS)
				report_MPI_ERROR(result);
		}
	}
#ifdef MPI_DEBUG
	std::cout << "relion_MPI_Allreduce: count = " << count << " datatype size = " << unitsize << " comm = " << comm << std::endl;
#endif

	return result;
}
#endif

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
		std::cout << " + Number of MPI processes                 = " << node.size << std::endl;
		if (nthreads > 1)
		{
        std::cout << " + Number of threads per MPI process       = " << nthreads << std::endl;
		std::cout << " + Total number of threads therefore       = " << nthreads * node.size << std::endl;
		}
		std::cout << " + Leader      (0) runs on host            = " << node.getHostName() << std::endl;
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
			std::cout << "  runs on host            = " << node.getHostName() << std::endl;
			std::cout.flush();
		}
		node.barrierWait();
	}

	if (node.isLeader())
	{
		std::cout << " ==========================" << std::endl;
	}
	std::cout.flush();

	// Try to flush all std::cout of all MPI processes before proceeding...
	sleep(1);
	node.barrierWait();
}
