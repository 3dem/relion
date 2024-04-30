#include "ctf_refinement_mpi.h"
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

#define TIMING 0

using namespace gravis;
using namespace aberration;


CtfRefinementProgramMpi::CtfRefinementProgramMpi(int argc, char *argv[])
		: CtfRefinementProgram(argc, argv)
{
	run_from_MPI = true;

	// Define a new MpiNode
	node = new MpiNode(argc, argv);
	rank = node->rank;
	nodeCount = node->size;

	// Don't put any output to screen for mpi followers
	verbosity = (node->isLeader()) ? 1 : 0;
}

void CtfRefinementProgramMpi::run()
{
	parseInput();

	if (verbosity > 0)
	{
		Log::beginSection("Initialising");
	}

	RefinementProgram::init();

	if (node->isLeader())
	{
		initTempDirectories();
	}

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize, particleSet.getTiltSeriesPixelSize(0));

	if (verbosity > 0)
	{
		Log::endSection();
	}

	std::vector<std::vector<int>> tomoIndices = ParticleSet::splitEvenly(particles, nodeCount);
    if(verbosity > 0)
    {
        Log::beginSection("Parallel tasks will be distributed as follows");
        for (int i = 0; i < nodeCount; i++)
        {
            Log::print(" Rank " + ZIO::itoa(i) + " will process " + ZIO::itoa(tomoIndices[i].size()) + " tomograms");
        }
        Log::print(" Progress below is only given for the process on Rank 0 ...");
        Log::endSection();
    }

	processTomograms(tomoIndices[rank], aberrationsCache, true);

	MPI_Barrier(MPI_COMM_WORLD);

	if (node->isLeader())
	{
		finalise();
	}
}
