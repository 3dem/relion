
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
#include <src/jaz/image/stack_helper.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/zio.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string particleFn, tomoSetFn, outPath;
	double minFreqPx, lowPass;
	bool oppositeHalf, predictCTF, applyDose, writeObs;
	int threads, tomoIndex, particleIndex, boxSize;

	//ReferenceMap reference;
	TomoReferenceMap referenceMap;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);

		parser.addSection("General options");

		particleFn = parser.getOption("--i", "Input particle set");
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		tomoIndex = textToInteger(parser.getOption("--ti", "Tomogram index", "0"));
		particleIndex = textToInteger(parser.getOption("--pi", "Particle index", "0"));
		boxSize = textToInteger(parser.getOption("--b", "Box size"));

		referenceMap.read(parser);

		minFreqPx = textToDouble(parser.getOption("--min_freq", "Min. image frequency [px]", "0"));

		predictCTF = parser.checkOption("--predict_CTF", "Modulate prediction by CTF");
		oppositeHalf = parser.checkOption("--opposite_half", "Make prediction from the opposite subset");
		applyDose = !parser.checkOption("--no_dose", "Do not apply dose weighting");
		writeObs = parser.checkOption("--write_obs", "Write out the observed images as well");

		lowPass = textToDouble(parser.getOption("--lowpass", "Blur by a Gaussian with this sigma [px]", "-1"));

		threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		outPath = parser.getOption("--o", "Output path");

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	if (outPath[outPath.length()-1] != '/')
	{
		outPath = outPath + "/";
	}

	std::string command = "mkdir -p " + outPath;
	int res = system(command.c_str());


	TomogramSet tomogramSet(tomoSetFn);

	ParticleSet dataSet(particleFn);
	std::vector<std::vector<ParticleIndex>> particles = dataSet.splitByTomogram(tomogramSet);

	AberrationsCache aberrationsCache(dataSet.optTable, boxSize, dataSet.getOriginalPixelSize(0));

	referenceMap.load(boxSize);

	const int s = referenceMap.getBoxSize();
	const int sh = s/2 + 1;


	{
		const int pc = particles[tomoIndex].size();

		if (pc <= particleIndex)
		{
			REPORT_ERROR_STR("Only " << pc << " particles found in tomogram " << tomoIndex);
		}

		const ParticleIndex part_id(particleIndex);

		Tomogram tomogram = tomogramSet.loadTomogram(tomoIndex, true);
		tomogram.validateParticleOptics({part_id}, dataSet);

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);

		const int fc = tomogram.frameCount;

		if (writeObs)
		{
			BufferedImage<float> obsStack(s,s,fc);
			BufferedImage<float> sliceRS(s,s);
			BufferedImage<tComplex<float>> sliceFS(sh,s);

			const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(part_id, fc, tomogram.optics.pixelSize);

			for (int f = 0; f < fc; f++)
			{
				d4Matrix projCut;

				TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram, traj[f],
					sliceFS, projCut, 1, true);

				Centering::shiftInSitu(sliceFS);
				FFT::inverseFourierTransform(sliceFS, sliceRS);

				obsStack.getSliceRef(f).copyFrom(sliceRS);
			}

			obsStack.write(
				outPath +  "tomo_" + ZIO::itoa(tomoIndex) + "_part_" + ZIO::itoa(particleIndex) + "_obs.mrc",
				tomogram.optics.pixelSize);
		}

		BufferedImage<float> predStack(s,s,fc);

		for (int f = 0; f < fc; f++)
		{
			CTF ctf = tomogram.getCtf(f, dataSet.getPosition(part_id));
			RawImage<float> doseSlice = doseWeights.getSliceRef(f);

			BufferedImage<fComplex> prediction = Prediction::predictModulated(
				part_id, dataSet, tomogram.projectionMatrices[f], s,
				ctf, tomogram.optics.pixelSize, aberrationsCache,
				referenceMap.image_FS,
				oppositeHalf? Prediction::OppositeHalf : Prediction::OwnHalf,
				predictCTF?   Prediction::AmplitudeModulated : Prediction::Unmodulated,
				applyDose? &doseSlice : 0,
				Prediction::CtfScaled);

			Centering::shiftInSitu(prediction);

			BufferedImage<float> predReal(s,s);
			FFT::inverseFourierTransform(prediction, predReal);

			predStack.getSliceRef(f).copyFrom(predReal);
		}

		predStack.write(
			outPath +  "tomo_" + ZIO::itoa(tomoIndex)
					+ "_part_" + ZIO::itoa(particleIndex) + "_pred.mrc",
			tomogram.optics.pixelSize);
	}

	return RELION_EXIT_SUCCESS;
}
