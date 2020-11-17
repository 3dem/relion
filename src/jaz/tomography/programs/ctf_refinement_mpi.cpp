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
	// Define a new MpiNode
	node = new MpiNode(argc, argv);
	rank = node->rank;
	nodeCount = node->size;

	// Don't put any output to screen for mpi slaves
	verb = (node->isMaster()) ? 1 : 0;

	if (nodeCount < 2)
	{
		REPORT_ERROR("SubtomoProgramMpi::read: this program needs to be run with at least two MPI processes!");
	}
}

void CtfRefinementProgramMpi::run()
{
	if (verb > 0)
	{
		Log::beginSection("Initialising");
	}

	RefinementProgram::init();

	if (node->isMaster())
	{
		initTempDirectories();
	}

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);

	if (verb > 0)
	{
		Log::endSection();
	}

	std::vector<std::vector<int>> tomoIndices = ParticleSet::splitEvenly(particles, nodeCount);

	processTomograms(tomoIndices[rank], aberrationsCache, verb, false);

	MPI_Barrier(MPI_COMM_WORLD);
	if (node->isMaster())
	{
		finalise();
	}
}
