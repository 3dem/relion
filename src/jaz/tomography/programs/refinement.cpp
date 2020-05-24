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
#include <src/jaz/tomography/data_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>
#include <iostream>
#include <src/time.h>

#define TIMING 0


using namespace gravis;

RefinementProgram::RefinementProgram(int argc, char *argv[])
:	argc(argc), argv(argv)
{}

void RefinementProgram::_readParams(IOParser &parser)
{
	int gen_section = parser.addSection("General refinement options");
	
	catFn = parser.getOption("--i", "Input particle set");
	tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
	boxSize = textToInteger(parser.getOption("--b", "Box size", "384"));
	
	ref1Fn = parser.getOption("--ref1", "Reference map, half 1");
	ref2Fn = parser.getOption("--ref2", "Reference map, half 2");
	maskFn = parser.getOption("--mask", "Reference mask", "");
	fscFn = parser.getOption("--fsc", "Star file containing the FSC of the reference", "");
	
	useFscThresh = !parser.checkOption("--fsc_act", "Use the actual FSC as the frq. weight");
	fscThreshWidth = textToDouble(parser.getOption("--fsc_thresh_width", "Width of the frq. weight flank", "5"));
			
	first_frame = textToInteger(parser.getOption("--f0", "First frame", "0"));
	last_frame = textToInteger(parser.getOption("--f1", "Last frame", "-1"));
	
	motFn = parser.getOption("--mot", "Particle trajectories", "");
	
	diag = parser.checkOption("--diag", "Write out diagnostic information");
	timing = parser.checkOption("--time", "Measure the elapsed time");
	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
	outDir = parser.getOption("--o", "Output filename pattern");
}

void RefinementProgram::init()
{
	if (outDir[outDir.length()-1] != '/')
	{
		outDir = outDir + "/";
	}
	
	int res = system(("mkdir -p "+outDir).c_str());
	
	{
		std::ofstream ofs(outDir+"/note.txt");
		
		ofs << "Command:\n\n";
		
		for (int i = 0; i < argc; i++)
		{
			ofs << argv[i] << ' ';
		}
		
		ofs << '\n';
	}

	tomogramSet = TomogramSet(tomoSetFn);
	
	dataSet = ParticleSet::load(catFn, motFn);
	particles = dataSet->splitByTomogram(tomogramSet);
		
	referenceMap = TomoReferenceMap(ref1Fn, ref2Fn, boxSize, maskFn, fscFn, useFscThresh, fscThreshWidth);
}

BufferedImage<float> RefinementProgram::computeFrequencyWeights(
		const Tomogram& tomogram,
		bool whiten, double sig2RampPower, double hiPass_px,
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
		frqWghts.getSliceRef(f) *= referenceMap.freqWeight;
	}
	
	BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);
	frqWghts *= doseWeights;
		
	return frqWghts;
}

