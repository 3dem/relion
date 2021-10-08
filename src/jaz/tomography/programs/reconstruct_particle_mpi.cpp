#include "reconstruct_particle_mpi.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/symmetry.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/time.h>
#include <iostream>


using namespace gravis;


void ReconstructParticleProgramMpi::readParameters(int argc, char *argv[])
{
	run_from_MPI = true;

	// Define a new MpiNode
	node = new MpiNode(argc, argv);
	rank = node->rank;
	nodeCount = node->size;

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	readBasicParameters(argc, argv);

	if (nodeCount < 2)
	{
		REPORT_ERROR("ReconstructParticleProgramMpi::read: this program needs to be run with at least two MPI processes!");
	}

	// Print out MPI info
	printMpiNodesMachineNames(*node);

	if (rank == 0)
	{
		outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);
	}
	else
	{
		outDir = ZIO::ensureEndingSlash(outDir);
	}

	ZIO::makeDir(outDir + "temp");
	tmpOutRootBase = outDir + "temp/sum_rank_";
	tmpOutRoot = tmpOutRootBase + ZIO::itoa(rank) + "_";
}

void ReconstructParticleProgramMpi::run()
{
	if (verb > 0)
	{
		Log::beginSection("Initialising");
	}

	TomogramSet tomoSet(optimisationSet.tomograms, verb > 0);
	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories, verb > 0);

	std::vector<std::vector<ParticleIndex>> particles = particleSet.splitByTomogram(tomoSet, verb > 0);

	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;

	const int s02D = (int)(binning * s + 0.5);

	const bool flip_value = true;
	const bool do_ctf = true;

	Tomogram tomo0 = tomoSet.loadTomogram(0, false);
	const double binnedOutPixelSize = tomo0.optics.pixelSize * binning;


	const long int voxelNum = (long int) sh * (long int) s * (long int) s;

	const double GB_per_thread =
			2.0 * voxelNum * 3.0 * sizeof(double)   // two halves  *  box size  *  (data (x2) + ctf)
			/ (1024.0 * 1024.0 * 1024.0);           // in GB

	if (max_mem_GB > 0)
	{
		const double maxThreads = max_mem_GB / GB_per_thread;

		if (maxThreads < outer_threads)
		{
			int lastOuterThreads = outer_threads;
			outer_threads = (int) maxThreads;

			if (verb > 0)
			{
				Log::print("Outer thread number reduced from " + ZIO::itoa(lastOuterThreads) +
					  " to " + ZIO::itoa(outer_threads) + " due to memory constraints (--mem).");
			}
		}
	}

	const int outCount = 2 * outer_threads; // One more count for the accumulated sum

	if (verb > 0)
	{
		Log::print("Memory required for accumulation: " + ZIO::itoa(GB_per_thread  * outer_threads) + " GB");
	}

	std::vector<BufferedImage<double>> ctfImgFS(outCount);
	std::vector<BufferedImage<dComplex>> dataImgFS(outCount);

	for (int i = 0; i < outCount; i++)
	{
		dataImgFS[i] = BufferedImage<dComplex>(sh,s,s);
		ctfImgFS[i] = BufferedImage<double>(sh,s,s),

		dataImgFS[i].fill(dComplex(0.0, 0.0));
		ctfImgFS[i].fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize, binnedOutPixelSize);

	if (verb > 0)
	{
		Log::endSection();
	}

	std::vector<std::vector<int>> tomoIndices = ParticleSet::splitEvenly(particles, nodeCount);

	processTomograms(
		tomoIndices[rank], tomoSet, particleSet, particles, aberrationsCache,
		dataImgFS, ctfImgFS, binnedOutPixelSize,
		s02D, do_ctf, flip_value, verb, false);

	std::vector<BufferedImage<double>> sumCtfImgFS(2);
	std::vector<BufferedImage<dComplex>> sumDataImgFS(2);

	size_t sizeData = sh*s*s;

	if (node->isLeader())
	{
		for (int i = 0; i < 2; i++)
		{
			sumDataImgFS[i] = BufferedImage<dComplex>(sh,s,s);
			sumCtfImgFS[i] = BufferedImage<double>(sh,s,s);
		}
	}

	for (int i = 0; i < 2; i++)
	{
		if (node->isLeader())
		{
			Log::print("Gathering half "+ZIO::itoa(i+1));
		}

		MPI_Reduce(dataImgFS[i].data, sumDataImgFS[i].data, sizeData,
				MY_MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce(ctfImgFS[i].data, sumCtfImgFS[i].data, sizeData,
				MY_MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		if (no_reconstruction) return;

		finalise(sumDataImgFS, sumCtfImgFS, binnedOutPixelSize);
	}

	// Delete temporary files (or try to; no error checking intentional)
	if (system(("rm -rf "+ tmpOutRootBase + "*.mrc").c_str()))
		std::cerr << "WARNING: deleting temporary files in folder " << tmpOutRootBase <<
		" failed." << std::endl;
}

