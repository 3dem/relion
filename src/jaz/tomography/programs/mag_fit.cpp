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
#include <src/jaz/tomography/projection/fwd_projection.h>

#include <src/jaz/math/Euler_angles_relion.h>
#include <src/euler.h>


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

	for (int t = 1; t < tc; t++)
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

		/*{
			const double pixelSize0 = tomogram.optics.pixelSize;

			for (int p = 0; p < 5; p++)
			{
				const ParticleIndex particle_id = particles[t][p];

				const d3Vector pos = dataSet.getPosition(particle_id);

				const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
							particle_id, fc, pixelSize0);


				const d4Matrix particleToTomo = dataSet.getMatrix4x4(
						particle_id, s, s, s);

				const int f = 19;

				BufferedImage<fComplex> observation(sh,s);

				d4Matrix projCut;

				TomoExtraction::extractFrameAt3D_Fourier(
						tomogram.stack, f, s, 1.0, tomogram.projectionMatrices[f],
						traj[f], observation, projCut, 1, false, true);

				std::vector<double> mags;

				for (double mag = 0.925; mag <= 1.075; mag += 0.025)
				{
					mags.push_back(mag);
				}

				BufferedImage<float> predictions(s,s,mags.size());
				BufferedImage<float> gradientsX(s,s,mags.size());
				BufferedImage<float> gradientsY(s,s,mags.size());

				const int hs = dataSet.getHalfSet(particle_id);

				for (int i = 0; i < mags.size(); i++)
				{
					const double mag = mags[i];

					d4Matrix projPart = mag * projCut * particleToTomo;

					BufferedImage<fComplex> prediction(sh,s);

					ForwardProjection::forwardProject(
							referenceMap.image_FS[hs], {projPart}, prediction, 1);

					BufferedImage<t2Vector<fComplex>> predGradient(sh,s);

					ForwardProjection::forwardProject2DGradient(
							referenceMap.image_FS[hs], {projPart}, predGradient, 1);

					BufferedImage<float> predRS;
					FFT::inverseFourierTransform(prediction, predRS);
					predictions.getSliceRef(i).copyFrom(predRS);


					BufferedImage<fComplex> gradFS_X(sh,s), gradFS_Y(sh,s);

					for (int y = 0; y < s;  y++)
					for (int x = 0; x < sh; x++)
					{
						gradFS_X(x,y) = predGradient(x,y).x;
						gradFS_Y(x,y) = predGradient(x,y).y;
					}

					BufferedImage<float> gradRS;
					FFT::inverseFourierTransform(gradFS_X, gradRS);
					gradientsX.getSliceRef(i).copyFrom(gradRS);

					FFT::inverseFourierTransform(gradFS_Y, gradRS);
					gradientsY.getSliceRef(i).copyFrom(gradRS);
				}

				predictions.write(outDir + "p_" + ZIO::itoa(particle_id.value) + "_prediction.mrc");
				gradientsX.write(outDir + "p_" + ZIO::itoa(particle_id.value) + "_gradientsX.mrc");
				gradientsY.write(outDir + "p_" + ZIO::itoa(particle_id.value) + "_gradientsY.mrc");

				BufferedImage<float> observation_RS;
				FFT::inverseFourierTransform(observation, observation_RS);

				observation_RS = ImageFilter::Gauss2D(observation_RS, 0, 3, true);
				observation_RS *= -1.f;

				observation_RS.write(outDir + "p_" + ZIO::itoa(particle_id.value) + "_observation.mrc");
			}

			std::exit(0);
		}*/

		/*TomoIsoMagFit isoFit(
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
		}*/




		/*TomoAnisoMagFit anisoFit(
			particles[t],
			tomogram,
			dataSet,
			referenceMap,
			freqWeights,
			boxSize,
			first_frame,
			last_frame,
			num_threads);



		BufferedImage<Equation2x2> equations = anisoFit.computeEquations();
		d2Matrix magMatrix = MagnificationHelper::solveLinearly(equations);

		std::cout << magMatrix << std::endl;*/




		d3Vector centre_of_mass = tomogram.computeCentreOfMass(
					dataSet, particles[t]);

		const double mean_depth = tomogram.getDepthOffset(
					tomogram.frameSequence[0], centre_of_mass);

		std::cout << "mean depth: " << mean_depth << std::endl;

		std::vector<ParticleIndex> particles_front, particles_back;

		for (int p = 0; p < pc; p++)
		{
			const ParticleIndex part_id = particles[t][p];

			const d3Vector pos = dataSet.getPosition(part_id);
			const double depth = tomogram.getDepthOffset(
					tomogram.frameSequence[0], pos);

			if (depth < mean_depth)
			{
				particles_front.push_back(part_id);
			}
			else
			{
				particles_back.push_back(part_id);
			}
		}

		std::cout << particles_front.size() << " in front\n";
		std::cout << particles_back.size() << " behind\n";

		TomoAnisoMagFit anisoFit_front(
			particles_front,
			tomogram,
			dataSet,
			referenceMap,
			freqWeights,
			boxSize,
			first_frame,
			last_frame,
			num_threads);

		TomoAnisoMagFit anisoFit_back(
			particles_back,
			tomogram,
			dataSet,
			referenceMap,
			freqWeights,
			boxSize,
			first_frame,
			last_frame,
			num_threads);

		BufferedImage<Equation2x2> equations_front = anisoFit_front.computeEquations();
		d2Matrix magMatrix_front = MagnificationHelper::solveLinearly(equations_front);

		std::cout << "front:\n" << magMatrix_front << std::endl;

		BufferedImage<Equation2x2> equations_back = anisoFit_back.computeEquations();
		d2Matrix magMatrix_back = MagnificationHelper::solveLinearly(equations_back);

		std::cout << "back:\n" << magMatrix_back << std::endl;



		/*BufferedImage<Equation2x2> equations = anisoFit.computeEquations();
		d2Matrix magMatrix = MagnificationHelper::solveLinearly(equations);

		std::cout << magMatrix << std::endl;*/



		/*d2Matrix identity;
		const double L2_0 = anisoFit.evaluateMag(identity);

		BufferedImage<Equation2x2> equations = anisoFit.computeEquations();

		d2Matrix magMatrix = MagnificationHelper::solveLinearly(equations);
		const double L2_1 = anisoFit.evaluateMag(magMatrix);

		std::cout << magMatrix << std::endl;
		std::cout << L2_0 << " -> " << L2_1 << " (" << (L2_0 - L2_1) << ")" << std::endl;*/



		/*for (double mx = 0.95; mx < 1.05; mx += 0.01)
		{
			const d2Matrix M(mx, 0, 0, 1);

			const double L2 = anisoFit.evaluateMag(M);

			std::cout << mx << ": " << L2 << std::endl;
		}*/



		/*std::vector<BufferedImage<Equation2x2>> equations = anisoFit.computeEquations_even_odd();

		std::vector<std::string> names {"all", "even", "odd"};

		for (int i = 0; i < 3; i++)
		{
			d2Matrix magMatrix = MagnificationHelper::solveLinearly(equations[i]);
			std::cout << names[i] << ":\n" << magMatrix << std::endl;
		}*/

		Log::endSection();
	}
}

