#include <src/args.h>
#include <mpi.h>
#include <src/ctf.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/motion_fit.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/motion/trajectory_set.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/index_sort.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/math/fcc.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/programs/align_mpi.h>

using namespace gravis;

AlignProgramMpi::AlignProgramMpi(int argc, char *argv[])
		: AlignProgram(argc, argv)
{
	run_from_MPI = true;

	// Define a new MpiNode
	node = new MpiNode(argc, argv);
	rank = node->rank;
	nodeCount = node->size;

	// Don't put any output to screen for mpi followers
	verbosity = (node->isLeader()) ? 1 : 0;
}

void AlignProgramMpi::run()
{
	parseInput();

	if (verbosity > 0)
	{
		Log::beginSection("Initialising");
	}

	initialise();

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize, particleSet.getOriginalPixelSize(0));

	if (verbosity > 0)
	{
		Log::endSection();
	}

	std::vector<std::vector<int>> tomoIndices = ParticleSet::splitEvenly(particles, nodeCount);

	processTomograms(tomoIndices[rank], aberrationsCache, false);

	MPI_Barrier(MPI_COMM_WORLD);

	if (node->isLeader())
	{
		finalise();
	}
}
