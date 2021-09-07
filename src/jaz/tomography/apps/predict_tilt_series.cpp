#include <src/args.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/util/log.h>
#include <src/args.h>

#include <omp.h>

using namespace gravis;


void run(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	try
	{
		run(argc, argv);
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

void run(int argc, char *argv[])
{
	OptimisationSet optimisationSet;

	IOParser parser;

	parser.setCommandLine(argc, argv);

	optimisationSet.read(
		parser,
		true,             // optimisation set
		true,    true,    // particles
		true,    true,    // tomograms
		true,	 false,   // trajectories
		false,   false,   // manifolds
		true,    true);   // reference

	const int gen_section = parser.addSection("General options");
	const int box_size = textToInteger(parser.getOption("--b", "Box size"));
	const std::string tomo_name = parser.getOption("--tn", "Tomogram name");
	const double spacing = textToDouble(parser.getOption("--bin", "Output binning", "1.0"));
	const bool do_subtract = parser.checkOption("--sub", "Subtract predicted tilt series from the original");
	const int n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
	const std::string outFn = parser.getOption("--o", "Output filename");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	ZIO::ensureParentDir(outFn);


	TomogramSet tomogram_set(optimisationSet.tomograms);
	const int tomo_index = tomogram_set.getTomogramIndex(tomo_name);

	if (tomo_index < 0)
	{
		REPORT_ERROR("There is no tomogram named " + tomo_name);
	}


	const bool load_image_data = do_subtract;

	Tomogram tomogram = tomogram_set.loadTomogram(tomo_index, load_image_data);

	const int w0 = tomogram.imageSize.x;
	const int h0 = tomogram.imageSize.y;
	const int fc = tomogram.frameCount;

	i2Vector binned_size = Resampling::getFourierCroppedSize2D(w0, h0, spacing, true);
	const int w1 = binned_size.x;
	const int h1 = binned_size.y;


	BufferedImage<float> new_stack(w1,h1,fc);

	ParticleSet particle_set(optimisationSet.particles, optimisationSet.trajectories, true);
	std::vector<std::vector<ParticleIndex>> all_particles = particle_set.splitByTomogram(tomogram_set, true);
	std::vector<ParticleIndex> particles = all_particles[tomo_index];


	tomogram.validateParticleOptics(particles, particle_set);

	AberrationsCache aberrations_cache(particle_set.optTable, box_size, particle_set.getOriginalPixelSize(0));
	BufferedImage<float> doseWeights = tomogram.computeDoseWeight(box_size, 1);


	TomoReferenceMap reference_map;
	reference_map.read(optimisationSet);
	reference_map.load(box_size, 1);


	Log::beginProgress("Predicting images", (int)(ceil(fc / (double)n_threads)));

	#pragma omp parallel for num_threads(n_threads)
	for (int f = 0; f < fc; f++)
	{
		const int th = omp_get_thread_num();

		if (th == 0)
		{
			Log::updateProgress(f);
		}

		BufferedImage<float> full_slice(w0,h0);

		const RawImage<float> doseSlice = doseWeights.getConstSliceRef(f);

		Prediction::predictMicrograph(
			f, particle_set, particles, tomogram,
			aberrations_cache, reference_map, &doseSlice,
			full_slice,
			Prediction::OwnHalf,
			Prediction::AmplitudeAndPhaseModulated,
			Prediction::CtfScaled);

		if (do_subtract)
		{
			full_slice -= tomogram.stack.getConstSliceRef(f);
			full_slice *= -1.f;
		}

		if (spacing != 1.0)
		{
			BufferedImage<float> binned_slice = Resampling::FourierCrop_fullStack(
						full_slice, spacing, n_threads, true);

			new_stack.getSliceRef(f).copyFrom(binned_slice);
		}
		else
		{
			new_stack.getSliceRef(f).copyFrom(full_slice);
		}
	}

	Log::endProgress();

	new_stack.write(outFn);
}
