#include "refinement.h"
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/optics/ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>
#include <iostream>
#include <src/time.h>

#define TIMING 0


using namespace gravis;

RefinementProgram::RefinementProgram(int argc, char *argv[])
:	argc(argc), argv(argv), run_from_MPI(false), verbosity(1)
{}

void RefinementProgram::_readParams(IOParser &parser)
{
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

	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Only process tomograms for which no temporary output data can be found");

	referenceMap.read(optimisationSet);

	static_noise = !parser.checkOption("--per_frame_noise", "Assume a different noise distribution for each frame");
	
	diag = parser.checkOption("--diag", "Write out diagnostic information");
	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
	outDir = parser.getOption("--o", "Output directory");

	run_from_GUI = is_under_pipeline_control();
}

void RefinementProgram::init()
{
	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);

	tomogramSet = TomogramSet(optimisationSet.tomograms, verbosity > 0);
	
	particleSet = ParticleSet(optimisationSet.particles, optimisationSet.trajectories, verbosity > 0);
	particles = particleSet.splitByTomogram(tomogramSet, verbosity > 0);
		
	referenceMap.load(boxSize, verbosity > 0);
}

BufferedImage<float> RefinementProgram::computeFrequencyWeights(
		const Tomogram& tomogram,
		bool whiten, double sig2RampPower, double hiPass_px, bool applyDoseWeight,
		int num_threads)
{
	const int s = boxSize;
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
				
				if (sig2RampPower != 0.0)
				{
					frqWghts1D[i] *= pow(i/(double)powSpec1D.size(), sig2RampPower);
				}
			}
			
			RawImage<float> fw = frqWghts.getSliceRef(f);
			RadialAvg::toFftwHalf_2D_lin(frqWghts1D, sh, s, fw);
		}

		if (static_noise)
		{
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
	}
	else
	{
		frqWghts.fill(1.0);
	}
	
	if (hiPass_px > 0.0)
	{
		const double b = -0.5 / (hiPass_px * hiPass_px);
		
		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < fc; f++)
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double xx = x;
			const double yy = y < s/2? y : y - s;
			const double r2 = xx*xx + yy*yy;
			
			frqWghts(x,y,f) *= 1.0 - exp(b*r2);
		}
	}

	for (int f = 0; f < fc; f++)
	{
		referenceMap.contributeWeight<float>(frqWghts.getSliceRef(f));
	}

	if (applyDoseWeight)
	{
		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);
		frqWghts *= doseWeights;
	}
		
	return frqWghts;
}

BufferedImage<int> RefinementProgram::findXRanges(
		const RawImage<float>& freqWeights,
		const RawImage<float>& doseWeights,
		double cutoffFraction)
{
	const int sh = freqWeights.xdim;
	const int s = freqWeights.ydim;
	const int fc = doseWeights.zdim;

	BufferedImage<int> out(s,fc);

	for (int f = 0; f < fc; f++)
	{
		double avg_freq = 0.0;

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			avg_freq += freqWeights(x,y,f);
		}

		avg_freq /= (s*sh);

		for (int y = 0; y < s; y++)
		{
			out(y,f) = 0;

			const double yy = y < s/2? y : y - s;
			const double xmax = sqrt(s*s/4 - yy*yy);

			for (int x = 0; x < sh && x <= xmax; x++)
			{
				const float fw = freqWeights(x,y,f) / avg_freq;
				const float dw = doseWeights(x,y,f);

				if (fw > cutoffFraction && dw > cutoffFraction)
				{
					out(y,f) = x+1;
				}
			}
		}
	}

	return out;
}

