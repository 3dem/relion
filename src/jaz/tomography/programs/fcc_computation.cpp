#include "fcc_computation.h"
#include <src/ctf.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/math/fcc.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <omp.h>


using namespace gravis;


FccProgram::FccProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
}

void FccProgram::readParams()
{
	IOParser parser;
	parser.setCommandLine(argc, argv);

	_readParams(parser);

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
}

void FccProgram::run()
{
	readParams();

	RefinementProgram::init();
	
	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;
	const bool flip_value = true;

	
	std::vector<double> pixelSizes(tc, 0.0);

	int fc_max = tomogramSet.getMaxFrameCount();

	BufferedImage<double> FCC3_sum(sh, fc_max, 3);
	FCC3_sum.fill(0.0);

	BufferedImage<double> FCC_by_t(sh, tc, 1);
	BufferedImage<double> all_FCCs(sh, fc_max, tc);

	std::ofstream ofs0(outDir + "all_FCCs.dat");

	
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		tomogram.validateParticleOptics(particles[t], particleSet);

		pixelSizes[t] = tomogram.optics.pixelSize;

		std::string tag = outDir + tomogram.name;


		BufferedImage<double> FCC3;
		const std::string tempFile = tag + "_FCC3.mrc";

		if (only_do_unfinished && ZIO::fileExists(tempFile))
		{
			FCC3.read(tag + "_FCC3.mrc");
		}
		else
		{
			FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				flip_value, num_threads);
			
			FCC3.write(tag + "_FCC3.mrc");
		}

		BufferedImage<double> fcc = FCC::divide(FCC3);
		fcc.write(tag + "_FCC.mrc");

		FCC3_sum += FCC3;

		BufferedImage<double> FCC3s = FCC::sumOverTime(FCC3);
		BufferedImage<double> FCCs = FCC::divide(FCC3s);

		const int fc = tomogram.frameCount;

		for (int x = 0; x < sh; x++)
		{
			FCC_by_t(x,t) = FCCs(x,0);

			for (int f = 0; f < fc; f++)
			{
				all_FCCs(x,f,t) = fcc(x,f);

				for (int d = 0; d < 3; d++)
				{
					FCC3_sum(x,f,d) += FCC3(x,f,d);
				}
			}
		}

		for (int x = 0; x < sh; x++)
		{
			ofs0 << x << ' ' << FCCs(x,0,0) << '\n';
		}

		ofs0 << '\n';

		{
			BufferedImage<double> scaleFactor(sh,fc);

			for (int f = 0; f < fc; f++)
			{
				scaleFactor(0,f) = 0.0;

				for (int x = 1; x < sh; x++)
				{
					scaleFactor(x,f) = FCC3(x,f,0) / FCC3(x,f,1);
				}
			}

			scaleFactor.write(tag + "_scaleFactor.mrc");
		}
	}

	ofs0.close();

	BufferedImage<double> FCC_sum = FCC::divide(FCC3_sum);

	FCC_sum.write(outDir + "FCC_by_frame.mrc");
	FCC_by_t.write(outDir + "FCC_by_tomo.mrc");
	all_FCCs.write(outDir + "FCC_by_frame_by_tomo.mrc");


	BufferedImage<double> FCC3_1D = FCC::divide(FCC::sumOverTime(FCC3_sum));

	std::ofstream ofs(outDir + "FCC.dat");

	for (int x = 0; x < sh; x++)
	{
		ofs << x << ' ' << FCC3_1D(x,0,0) << '\n';
	}

	ofs.close();
}

