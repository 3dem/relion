#include "mag_fit.h"
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/math/Zernike_helper.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/optics/tomo_mag_fit.h>
#include <omp.h>

using namespace gravis;
using namespace aberration;


MagFitProgram::MagFitProgram(int argc, char *argv[])
:	RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void MagFitProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);

		int defocus_section = parser.addSection("Alignment options");


		initial_step = textToDouble(parser.getOption("--ins", "Initial step (in %)", "2"));


		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			parser.writeUsage(std::cout);
			std::exit(-1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void MagFitProgram::run()
{
	Log::beginSection("Initialising");

	RefinementProgram::init();

	const int s = boxSize;
	const int sh = s/2 + 1;
	const int tc = particles.size();
	const int gc = dataSet.numberOfOpticsGroups();
	const bool flip_value = true;

	Log::endSection();

	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));
		Log::print("Loading");

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);

		const int fc = tomogram.frameCount;

		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, num_threads);

		dataSet.checkTrajectoryLengths(
				particles[t][0], pc, fc, "MagFitProgram::run");

		TomoIsoMagFit isoFit(
			particles[t],
			tomogram,
			dataSet,
			referenceMap,
			freqWeights,
			boxSize,
			first_frame,
			last_frame,
			num_threads);


		std::vector<d3Vector> results_full(0);
		std::vector<d3Vector> results_defocus_only(0);
		std::vector<d3Vector> results_scale_only(0);

		for (double mag = 0.925; mag <= 1.075; mag += 0.025)
		{
			std::cout << mag << std::endl;

			d2Vector full = isoFit.computeErrorAndSlope(mag, true, true);
			d2Vector defocus_only = isoFit.computeErrorAndSlope(mag, false, true);
			d2Vector scale_only = isoFit.computeErrorAndSlope(mag, true, false);

			results_full.push_back(d3Vector(mag, full[0], full[1]));
			results_defocus_only.push_back(d3Vector(mag, defocus_only[0], defocus_only[1]));
			results_scale_only.push_back(d3Vector(mag, scale_only[0], scale_only[1]));
		}

		std::ofstream file_full(outDir + "t_" + ZIO::itoa(t) + "_full.dat");
		std::ofstream file_defocus_only(outDir + "t_" + ZIO::itoa(t) + "_defocus_only.dat");
		std::ofstream file_scale_only(outDir + "t_" + ZIO::itoa(t) + "_scale_only.dat");

		file_full.precision(16);
		file_defocus_only.precision(16);
		file_scale_only.precision(16);

		for (int i = 0; i < results_full.size(); i++)
		{
			file_full << results_full[i][0] << ' ' << results_full[i][1] << '\n';
			file_defocus_only << results_defocus_only[i][0] << ' ' << results_defocus_only[i][1] << '\n';
			file_scale_only << results_scale_only[i][0] << ' ' << results_scale_only[i][1] << '\n';
		}

		file_full << '\n';
		file_defocus_only << '\n';
		file_scale_only << '\n';

		for (int i = 0; i < results_full.size(); i++)
		{
			file_full << results_full[i][0] << ' ' << results_full[i][2] << '\n';
			file_defocus_only << results_defocus_only[i][0] << ' ' << results_defocus_only[i][2] << '\n';
			file_scale_only << results_scale_only[i][0] << ' ' << results_scale_only[i][2] << '\n';
		}

		Log::endSection();
	}
}

