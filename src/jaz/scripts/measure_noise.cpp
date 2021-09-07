
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/new_ft.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/reference_map.h>
#include <src/jaz/single_particle/new_reference_map.h>
#include <src/jaz/tomography/reference_map.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/util/log.h>
#include <src/jaz/image/stack_helper.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/zio.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string outDir;
	int num_threads, tomoIndex, boxSize;

	OptimisationSet optimisationSet;

	//ReferenceMap reference;
	TomoReferenceMap referenceMap;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);

		optimisationSet.read(
			parser,
			true,           // optimisation set
			true,   true,   // particles
			true,   true,   // tomograms
			true,   false,  // trajectories
			false,  false,  // manifolds
			true,   true);  // reference

		int gen_section = parser.addSection("General refinement options");

		boxSize = textToInteger(parser.getOption("--b", "Box size"));

		referenceMap.read(optimisationSet);

		tomoIndex = textToInteger(parser.getOption("--ti", "Tomogram index", "0"));
		boxSize = textToInteger(parser.getOption("--b", "Box size"));

		num_threads = textToInteger(parser.getOption("--j", "Number of threads", "8"));
		outDir = parser.getOption("--o", "Output path");

		parser.checkForErrors();

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
		return RELION_EXIT_FAILURE;
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);


	TomogramSet tomogramSet(optimisationSet.tomograms);

	ParticleSet dataSet(optimisationSet.particles);
	std::vector<std::vector<ParticleIndex>> particles = dataSet.splitByTomogram(tomogramSet);

	AberrationsCache aberrationsCache(dataSet.optTable, boxSize, dataSet.getOriginalPixelSize(0));

	referenceMap.load(boxSize);


	Tomogram tomogram = tomogramSet.loadTomogram(tomoIndex, true);
	tomogram.validateParticleOptics(particles[tomoIndex], dataSet);

	BufferedImage<float> doseWeights = tomogram.computeDoseWeight(boxSize, 1.0);

	const int s = referenceMap.getBoxSize();
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;

	std::vector<std::vector<double>> sum_delta2_f(fc, std::vector<double>(sh, 0.0));
	std::vector<std::vector<double>> sum_obsPwr_f(fc, std::vector<double>(sh, 0.0));
	std::vector<std::vector<double>> sum_weight_f(fc, std::vector<double>(sh, 0.0));


	const int pc = particles[tomoIndex].size();

	Log::beginProgress("Measuring noise", pc);

	for (int p = 0; p < pc; p++)
	{
		Log::updateProgress(p);

		const ParticleIndex part_id(p);

		const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(part_id, fc, tomogram.optics.pixelSize);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < fc; f++)
		{
			d4Matrix projCut;

			BufferedImage<tComplex<float>> observation(sh,s);

			TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram, traj[f],
				observation, projCut, 1, true);

			CTF ctf = tomogram.getCtf(f, dataSet.getPosition(part_id));
			RawImage<float> doseSlice = doseWeights.getSliceRef(f);

			BufferedImage<fComplex> prediction = Prediction::predictModulated(
				part_id, dataSet, tomogram.projectionMatrices[f], s,
				ctf, tomogram.optics.pixelSize, aberrationsCache,
				referenceMap.image_FS,
				Prediction::OwnHalf,
				Prediction::AmplitudeModulated,
				&doseSlice,
				Prediction::CtfScaled);

			for (int y = 0; y < sh; y++)
			for (int x = 0; x < s;  x++)
			{
				const double xx = x;
				const double yy = y < s/2? y : y - s;
				const double r = sqrt(xx*xx + yy*yy);

				const fComplex delta = prediction(x,y) - (-observation(x,y));

				const int ri = (int) r;

				if (ri < sh)
				{
					sum_delta2_f[f][ri] += delta.norm();
					sum_obsPwr_f[f][ri] += observation(x,y).norm();
					sum_weight_f[f][ri] += 1.0;
				}
			}
		}
	}

	Log::endProgress();


	std::vector<double> sum_delta2(sh, 0.0);
	std::vector<double> sum_obsPwr(sh, 0.0);
	std::vector<double> sum_weight(sh, 0.0);

	for (int f = 0; f < fc; f++)
	{
		std::vector<double> var_delta(sh, 0.0), var_obs(sh, 0.0);

		for (int ri = 0; ri < sh; ri++)
		{
			if (sum_weight_f[f][ri] > 0.0)
			{
				var_delta[ri] = sum_delta2_f[f][ri] / sum_weight_f[f][ri];
				var_obs[ri] = sum_obsPwr_f[f][ri] / sum_weight_f[f][ri];
			}

			sum_delta2[ri] += sum_delta2_f[f][ri];
			sum_obsPwr[ri] += sum_obsPwr_f[f][ri];
			sum_weight[ri] += sum_weight_f[f][ri];
		}

		std::ofstream out(outDir + "noise_frame_" + ZIO::itoa(f) + ".dat");

		for (int ri = 0; ri < sh; ri++)
		{
			out << ri << ' ' << var_delta[ri] << '\n';
		}

		out << '\n';

		for (int ri = 0; ri < sh; ri++)
		{
			out << ri << ' ' << var_obs[ri] << '\n';
		}


	}

	std::vector<double> var_delta(sh, 0.0), var_obs(sh, 0.0);

	for (int ri = 0; ri < sh; ri++)
	{
		if (sum_weight[ri] > 0.0)
		{
			var_delta[ri] = sum_delta2[ri] / sum_weight[ri];
			var_obs[ri] = sum_obsPwr[ri] / sum_weight[ri];
		}
	}

	std::ofstream out(outDir + "noise_all_frames.dat");

	for (int ri = 0; ri < sh; ri++)
	{
		out << ri << ' ' << var_delta[ri] << '\n';
	}

	out << '\n';

	for (int ri = 0; ri < sh; ri++)
	{
		out << ri << ' ' << var_obs[ri] << '\n';
	}

	return RELION_EXIT_SUCCESS;
}
