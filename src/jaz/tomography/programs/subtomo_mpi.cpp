#include "subtomo_mpi.h"
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/projection/point_insertion.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/time.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>

using namespace gravis;

void SubtomoProgramMpi::readParameters(int argc, char *argv[])
{
	run_from_MPI = true;

	// Define a new MpiNode
	node = new MpiNode(argc, argv);
	rank = node->rank;
	nodeCount = node->size;

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	IOParser parser;

	parser.setCommandLine(argc, argv);

	readBasicParameters(parser);
	do_sum_all = false;

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	if (nodeCount < 2)
	{
		REPORT_ERROR("SubtomoProgramMpi::read: this program needs to be run with at least two MPI processes!");
	}

	// Print out MPI info
	printMpiNodesMachineNames(*node);

	if (do_gridding_precorrection)
	{
		do_narrow_circle_crop = true;
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);

	if (node->isLeader())
	{
		ZIO::makeDir(outDir + "/Subtomograms");
	}
}

void SubtomoProgramMpi::run()
{
	TomogramSet tomogramSet(optimisationSet.tomograms, verb > 0);

	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories, verb > 0);
	std::vector<std::vector<ParticleIndex> > particles = particleSet.splitByTomogram(tomogramSet, verb > 0);

	if (cropSize < 0) cropSize = boxSize;

	bool do_ctf = true;

	const long int tc = particles.size();
	const long int s2D = boxSize;

	const long int s3D = cropSize;
	const long int sh3D = s3D / 2 + 1;

	const long int s02D = (int)(binning * s2D + 0.5);

	const double relative_box_scale = cropSize / (double) boxSize;
	const double binned_pixel_size = binning * particleSet.getOriginalPixelSize(0);

	if (node->isLeader())
	{
		initialise(particleSet, particles, tomogramSet);
	}

	MPI_Barrier(MPI_COMM_WORLD);


	BufferedImage<float> dummy;

	AberrationsCache aberrationsCache(particleSet.optTable, s2D, binned_pixel_size);


	std::vector<std::vector<int>> tomoIndices = ParticleSet::splitEvenly(particles, nodeCount);

	processTomograms(
			tomoIndices[rank],
			tomogramSet,
			particleSet,
			particles,
			aberrationsCache,
			s02D,
			s2D,
			s3D,
			relative_box_scale,
			do_ctf,
			verb,
			dummy,
			dummy);
}
