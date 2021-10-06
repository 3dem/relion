
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
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/stack_helper.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/zio.h>

#include <omp.h>

using namespace gravis;


BufferedImage<float> computeFrequencyWeights(
		const Tomogram& tomogram, TomoReferenceMap& referenceMap,
		bool whiten, int num_threads)
{
	const int s = referenceMap.image_FS[0].ydim;
	const int sh = s / 2 + 1;
	const int fc = tomogram.frameCount;

	BufferedImage<float> frqWghts(sh,s,fc);

	if (whiten)
	{
		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < fc; f++)
		{
			BufferedImage<double> powSpec = PowerSpectrum::periodogramAverage2D(
				tomogram.stack, s, s, 2.0, f, false);

			std::vector<double> powSpec1D = RadialAvg::fftwHalf_2D_lin(powSpec);

			std::vector<float> frqWghts1D(powSpec1D.size());

			for (int i = 0; i < powSpec1D.size(); i++)
			{
				frqWghts1D[i] = (float)(1.0 / powSpec1D[i]);
			}

			RawImage<float> fw = frqWghts.getSliceRef(f);
			RadialAvg::toFftwHalf_2D_lin(frqWghts1D, sh, s, fw);
		}


		BufferedImage<float> staticSigma2(sh,s);
		staticSigma2.fill(0.f);

		for (int f = 0; f < fc; f++)
		{
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const double val = frqWghts(x,y,f);

				if (val > 0.f)
				{
					staticSigma2(x,y) += 1.f / val;
				}
			}
		}

		for (int f = 0; f < fc; f++)
		{
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				frqWghts(x,y,f) = fc / staticSigma2(x,y);
			}
		}
	}
	else
	{
		frqWghts.fill(1.0);
	}

	for (int f = 0; f < fc; f++)
	{
		referenceMap.contributeWeight<float>(frqWghts.getSliceRef(f));
	}

	BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);
	frqWghts *= doseWeights;

	return frqWghts;
}



int main(int argc, char *argv[])
{
	std::string outDir;
	int num_threads, boxSize;

	OptimisationSet optimisationSet;

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


	const int f0 = 10;
	const int f1 = 30;

	for (int tomoIndex = 0; tomoIndex < tomogramSet.size(); tomoIndex++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(tomoIndex, true);
		tomogram.validateParticleOptics(particles[tomoIndex], dataSet);

		std::string outDirTomo = outDir + tomogram.name + "/";
		ZIO::ensureParentDir(outDirTomo);

		const int s = referenceMap.getBoxSize();
		const int sh = s/2 + 1;
		const int fc = tomogram.frameCount;

		std::vector<std::vector<double>> sum_prdObs_f(fc, std::vector<double>(sh, 0.0));
		std::vector<std::vector<double>> sum_prdSqr_f(fc, std::vector<double>(sh, 0.0));


		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);

		BufferedImage<float> frqWeight = computeFrequencyWeights(
			tomogram, referenceMap, true, num_threads);


		const int pc = particles[tomoIndex].size();

		Log::beginProgress("Measuring scale", pc);

		for (int p = 0; p < pc; p++)
		{
			Log::updateProgress(p);

			const ParticleIndex part_id = particles[tomoIndex][p];

			const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(part_id, fc, tomogram.optics.pixelSize);

			#pragma omp parallel for num_threads(num_threads)
			for (int f = f0; f <= f1; f++)
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
					Prediction::CtfUnscaled);


				const double cos_f = tomogram.projectionMatrices[f](2,2);
				const double sin2_f = 1.0 - cos_f * cos_f;
				const double tan2_f = sin2_f / (cos_f * cos_f);

				for (int y = 0; y < sh; y++)
				for (int x = 0; x < s;  x++)
				{
					const double xx = x;
					const double yy = y < s/2? y : y - s;
					const double r = sqrt(xx*xx + yy*yy);

					const fComplex obs = -observation(x,y);
					const fComplex prd =  prediction(x,y);

					const int ri = (int) r;

					if (ri < sh)
					{
						sum_prdObs_f[f][ri] += frqWeight(x,y) * (prd.real * obs.real + prd.imag * obs.imag);
						sum_prdSqr_f[f][ri] += frqWeight(x,y) * (prd.real * prd.real + prd.imag * prd.imag);
					}
				}
			}
		}

		Log::endProgress();


		std::vector<double> sum_prdObs_byFreq(sh, 0.0);
		std::vector<double> sum_prdSqr_byFreq(sh, 0.0);

		std::vector<double> sum_prdObs_byFrame(fc, 0.0);
		std::vector<double> sum_prdSqr_byFrame(fc, 0.0);
		std::vector<double> ratio_byFrame(fc, 0.0);

		double sum_prdObs_tomo = 0.0;
		double sum_prdSqr_tomo = 0.0;

		const double eps = 0.01;

		for (int f = 0; f < fc; f++)
		{
			std::vector<double> ratio(sh, 0.0);

			for (int ri = 0; ri < sh; ri++)
			{
				if (sum_prdSqr_f[f][ri] > eps)
				{
					ratio[ri] = sum_prdObs_f[f][ri] / sum_prdSqr_f[f][ri];
				}

				sum_prdObs_byFreq[ri] += sum_prdObs_f[f][ri];
				sum_prdSqr_byFreq[ri] += sum_prdSqr_f[f][ri];

				sum_prdObs_byFrame[f] += sum_prdObs_f[f][ri];
				sum_prdSqr_byFrame[f] += sum_prdSqr_f[f][ri];

				sum_prdObs_tomo += sum_prdObs_f[f][ri];
				sum_prdSqr_tomo += sum_prdSqr_f[f][ri];
			}

			std::ofstream out(outDirTomo + "scale_over_freq__frame_" + ZIO::itoa(f) + ".dat");

			for (int ri = 0; ri < sh; ri++)
			{
				if (ratio[ri] > 0.0)
				{
					out << ri << ' ' << ratio[ri] << '\n';
				}
			}

			if (sum_prdSqr_byFrame[f] > eps)
			{
				ratio_byFrame[f] = sum_prdObs_byFrame[f] / sum_prdSqr_byFrame[f];
			}
		}

		std::cout << tomogram.name << ": " << (sum_prdObs_tomo / sum_prdSqr_tomo) << '\n';

		std::vector<double> ratio_byFreq(sh, 0.0);
		std::ofstream out(outDirTomo + "scale_over_freq__all_frames.dat");

		for (int ri = 0; ri < sh; ri++)
		{
			if (sum_prdSqr_byFreq[ri] > eps)
			{
				ratio_byFreq[ri] = sum_prdObs_byFreq[ri] / sum_prdSqr_byFreq[ri];
				out << ri << ' ' << ratio_byFreq[ri] << '\n';
			}
		}


		std::ofstream out2(outDirTomo + "scale_over_frame.dat");
		std::ofstream out3(outDirTomo + "thickness_by_frame.dat");

		for (int f = 0; f < fc; f++)
		{
			if (ratio_byFrame[f] > 0.0)
			{
				out2 << f << ' ' << ratio_byFrame[f] << '\n';
			}

			const double cos_f = tomogram.projectionMatrices[f](2,2);
			const double sin2_f = 1.0 - cos_f * cos_f;
			const double tan2_f = sin2_f / (cos_f * cos_f);

			out3 << f << ' ' << 1.0 / sqrt(1 + tan2_f) << '\n';
		}
	}

	return RELION_EXIT_SUCCESS;
}
