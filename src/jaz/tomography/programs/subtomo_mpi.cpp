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
	// Define a new MpiNode
	node = new MpiNode(argc, argv);
	rank = node->rank;
	nodeCount = node->size;

	// Don't put any output to screen for mpi slaves
	verb = (node->isMaster()) ? verb : 0;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);

		readBasicParameters(parser);

		Log::readParams(parser);

		if (parser.checkForErrors()) std::exit(-1);

	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	if (nodeCount < 2)
		REPORT_ERROR("SubtomoProgramMpi::read ERROR: this program needs to be run with at least two MPI processes!");
	// Print out MPI info
	printMpiNodesMachineNames(*node);

	if (do_gridding_precorrection)
	{
		do_narrow_circle_crop = true;
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);

	if (node->isMaster())
	{
		int res = system(("mkdir -p " + outDir + "/Subtomograms").c_str());
	}
}

void SubtomoProgramMpi::run()
{
	TomogramSet tomogramSet(optimisationSet.tomograms);

	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories);
	std::vector<std::vector<ParticleIndex> > particles = particleSet.splitByTomogram(tomogramSet);

	if (cropSize < 0) cropSize = boxSize;

	bool do_ctf = true;

	const long int tc = particles.size();
	const long int s2D = boxSize;

	const long int s3D = cropSize;
	const long int sh3D = s3D / 2 + 1;

	const long int s02D = (int)(binning * s2D + 0.5);

	const double relative_box_scale = cropSize / (double) boxSize;

	if (node->isMaster())
	{
		writeParticleSet(particleSet, particles);
	}

	BufferedImage<float> sum_data, sum_weights;

	if (do_sum_all)
	{
		sum_data.resize(s3D,s3D,s3D);
		sum_data.fill(0.0);

		sum_weights.resize(sh3D,s3D,s3D);
		sum_weights.fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, s2D);

	// determine tomogram range based on node rank:
	const int first_tomo = rank * tc / nodeCount;
	const int last_tomo = (rank == nodeCount - 1)? tc - 1 : (rank + 1) * tc / nodeCount - 1;

	processTomograms(
			first_tomo,
			last_tomo,
			tomogramSet,
			particleSet,
			particles,
			aberrationsCache,
			s02D,
			s2D,
			s3D,
			relative_box_scale,
			do_ctf,
			1,
			sum_data,
			sum_weights);


	if (do_sum_all)
	{
		BufferedImage<float> sum_data_global, sum_weights_global;
		if (node->isMaster())
		{
			sum_data_global.resize(sum_data);
			sum_weights_global.resize(sum_weights);
		}
		MPI_Reduce(MULTIDIM_ARRAY(sum_data), MULTIDIM_ARRAY(sum_data_global), sum_data.getSize(),
					  MY_MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MULTIDIM_ARRAY(sum_weights), MULTIDIM_ARRAY(sum_weights_global), sum_weights.getSize(),
					  MY_MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (node->isMaster())
		{
			sum_data_global.write(outDir + "sum_data.mrc");
			sum_weights_global.write(outDir + "sum_weight.mrc");

			BufferedImage<float> dataImgDivRS(s3D, s3D, s3D);
			dataImgDivRS.fill(0.0);

			if (SNR > 0.0)
			{
				Reconstruction::ctfCorrect3D_Wiener(
						sum_data_global, sum_weights_global, dataImgDivRS,
						1.0 / SNR, num_threads);
			}
			else
			{
				Reconstruction::ctfCorrect3D_heuristic(
						sum_data_global, sum_weights_global, dataImgDivRS,
						0.001, num_threads);
			}

			dataImgDivRS.write(outDir + "sum_div.mrc");
		}
	}
}