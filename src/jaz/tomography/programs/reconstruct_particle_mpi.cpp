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


ReconstructParticleProgramMpi::ReconstructParticleProgramMpi(int rank, int nodeCount)
: rank(rank),
  nodeCount(nodeCount)
{
}

void ReconstructParticleProgramMpi::readParameters(int argc, char *argv[])
{
	readBasicParameters(argc, argv);

	if (rank == 0)
	{
		if (outTag.find_last_of("/") != std::string::npos)
		{
			std::string dir = outTag.substr(0, outTag.find_last_of("/"));
			int res = system(("mkdir -p "+dir).c_str());
		}

		{
			std::ofstream ofs(outTag+"_note.txt");

			ofs << "Command:\n\n";

			for (int i = 0; i < argc; i++)
			{
				ofs << argv[i] << ' ';
			}

			ofs << '\n';
		}
	}
}

void ReconstructParticleProgramMpi::run()
{
	Log::beginSection("Initialising");

	TomogramSet tomoSet(optimisationSet.tomograms);
	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories);

	std::vector<std::vector<ParticleIndex>> particles = particleSet.splitByTomogram(tomoSet);

	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;

	const int s02D = (int)(binning * s + 0.5);

	const bool flip_value = true;
	const bool do_ctf = true;

	Tomogram tomo0 = tomoSet.loadTomogram(0, false);
	const double binnedOutPixelSize = tomo0.optics.pixelSize * binning;


	const long int voxelNum = (long int) sh * (long int) s * (long int) s;

	if (max_mem_GB > 0)
	{
		const double GB_per_thread =
		   2.0 * voxelNum * 3.0 * sizeof(double)   // two halves  *  box size  *  (data (x2) + ctf)
			/ (1024.0 * 1024.0 * 1024.0);          // in GB

		const double maxThreads = max_mem_GB / GB_per_thread;

		if (maxThreads < outer_threads)
		{
			int lastOuterThreads = outer_threads;
			outer_threads = (int) maxThreads;

			Log::print("Outer thread number reduced from " + ZIO::itoa(lastOuterThreads) +
					  " to " + ZIO::itoa(outer_threads) + " due to memory constraints (--mem).");
		}
	}

	const int outCount = 2 * outer_threads;

	Log::print("Memory required for accumulation: " + ZIO::itoa(
			(3.0 * sizeof(double) * (long int) outCount * (double)voxelNum)
			  / (1024.0 * 1024.0 * 1024.0)
			) + " GB");

	std::vector<BufferedImage<double>> ctfImgFS(outCount), psfImgFS(outCount);
	std::vector<BufferedImage<dComplex>> dataImgFS(outCount);

	for (int i = 0; i < outCount; i++)
	{
		dataImgFS[i] = BufferedImage<dComplex>(sh,s,s);
		ctfImgFS[i] = BufferedImage<double>(sh,s,s),
		psfImgFS[i] = BufferedImage<double>(sh,s,s);

		dataImgFS[i].fill(dComplex(0.0, 0.0));
		ctfImgFS[i].fill(0.0);
		psfImgFS[i].fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);

	Log::endSection();

	// determine tomogram range based on node rank:
	const int first_tomo = rank * tc / nodeCount;
	const int last_tomo = (rank == nodeCount - 1)? tc - 1 : (rank + 1) * tc / nodeCount - 1;

	processTomograms(
		first_tomo, last_tomo, tomoSet, particleSet, particles, aberrationsCache,
		dataImgFS, ctfImgFS, psfImgFS, binnedOutPixelSize,
		s02D, do_ctf, flip_value, 1);


	if (outCount > 2)
	{
		Log::print("Merging volumes");

		for (int i = 2; i < outCount; i++)
		{
			dataImgFS[i%2] += dataImgFS[i];
			ctfImgFS[i%2] += ctfImgFS[i];

			if (explicit_gridding)
			{
				psfImgFS[i%2] += psfImgFS[i];
			}
		}
	}


	/*

	if (rank > 0)
	{
		write dataImgFS, ctfImgFS and psfImgFS to disk (may need to IFFT them first)
	}
	else
	{
		barrier
		read all results from disk and add them up
	}

	*/


	if (rank == 0)
	{
		if (no_reconstruction) return;

		if (symmName != "C1")
		{
			Log::print("Applying symmetry");

			for (int half = 0; half < 2; half++)
			{
				dataImgFS[half] = Symmetry::symmetrise_FS_complex(
							dataImgFS[half], symmName, num_threads);

				ctfImgFS[half] = Symmetry::symmetrise_FS_real(
							ctfImgFS[half], symmName, num_threads);

				if (explicit_gridding)
				{
					psfImgFS[half] = Symmetry::symmetrise_FS_real(
							psfImgFS[half], symmName, num_threads);
				}
			}
		}

		finalise(dataImgFS, ctfImgFS, psfImgFS, binnedOutPixelSize);
	}
}

