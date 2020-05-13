#include "bfactor_fit.h"
#include <src/ctf.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/math/fcc.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/tomography/data_set.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <omp.h>


using namespace gravis;


BfactorFitProgram::BfactorFitProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void BfactorFitProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);
		
		int def_section = parser.addSection("B factor refinement options");
		
		useL2 = parser.checkOption("--L2", "Use L2 based formulation instead of the FCC");
		kMin = textToInteger(parser.getOption("--kmin", "Min. freq. used for estimation [Px]", "20"));
		kMin2 = textToInteger(parser.getOption("--kmin2", "Throw out frames with an estimated sigma below this [Px]", "10"));
		useCache = parser.checkOption("--cache", "Use the FCC3 files computed during a previous run (if available)");
		
		if (parser.checkForErrors()) std::exit(-1);
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void BfactorFitProgram::run()
{
	REPORT_ERROR("B-factor fit is currently not operational. Reason: file format refactoring in progress.");
	RefinementProgram::init();
	
	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;
	const bool flip_value = true;
	
	
	TomogramSet tomogramSetOut = tomogramSet;
	
	std::vector<std::vector<double>> B_t(tc), k_t(tc), cumulativeDose(tc);
	
	std::vector<double> pixelSizes(tc, 0.0);
			
	
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		
		std::string tag = outDir + ZIO::itoa(t);
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		pixelSizes[t] = tomogram.optics.pixelSize;
		
		BufferedImage<double> FCC3;
		
		const std::string fcc3fn = tag + "_FCC3.mrc";
		
		if (useCache && exists(fcc3fn))
		{
			std::cout << "        reading: " << fcc3fn << std::endl;
			FCC3.read(fcc3fn);
		}
		else
		{
			FCC3 = FCC::compute3(
				dataSet, particles[t], tomogram, referenceMap.image_FS,
				flip_value, num_threads);
			
			FCC3.write(fcc3fn);
		}
		
		FCC::divide(FCC3).write(tag + "_FCC.mrc");
		
		BufferedImage<double> plot;
		
		std::pair<std::vector<double>, std::vector<double>> bkFacs = 
				Damage::fitBkFactors(FCC3, boxSize, tomogram.optics.pixelSize, kMin, sh, plot, useL2);
		
		plot.write(tag + "_fit.mrc");
		
		B_t[t] = bkFacs.first;
		k_t[t] = bkFacs.second;
		
		cumulativeDose[t] = tomogram.cumulativeDose;
				
	}
	
	Damage::renormalise(B_t, k_t);
	
	std::ofstream bfacsByDoseDat(outDir + "B_factors_by_dose.dat");
	std::ofstream kfacsByDoseDat(outDir + "scale_factors_by_dose.dat");
	
	
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		
		std::string tag = outDir + ZIO::itoa(t);
		
		BufferedImage<double> FCC3;		
		FCC3.read(tag + "_FCC3.mrc");
		
		const int fc = FCC3.ydim;
		
		
		BufferedImage<double> fit0 = Damage::renderFit(
					B_t[t], k_t[t], boxSize, pixelSizes[t], sh);
		
		fit0.write(tag + "_weights.mrc");
				
		
		/*MetaDataTable mdt;
		
		std::ofstream bfacsDat(tag+"_B_factors.dat");
		std::ofstream kfacsDat(tag+"_scale_factors.dat");
		
		for (int f = 0; f < fc; f++)
		{
			double B = B_t[t][f];
			double a = k_t[t][f];
						
			double k = a > 0.0? log(a) : -100.0;
			
			mdt.addObject();
			mdt.setValue(EMDL_MICROGRAPH_FRAME_NUMBER, f);
			mdt.setValue(EMDL_MICROGRAPH_PRE_EXPOSURE, cumulativeDose[t][f]);
			mdt.setValue(EMDL_CTF_BFACTOR, B);
			mdt.setValue(EMDL_CTF_SCALEFACTOR, a);
			
			if (a > 1e-3)
			{
				bfacsDat << f << ' ' << B << '\n';
				kfacsDat << f << ' ' << k << '\n';
				
				double dose = cumulativeDose[t][f];
				
				bfacsByDoseDat << dose << ' ' << B << '\n';
				kfacsByDoseDat << dose << ' ' << k << '\n';
			}
		}
		
		bfacsByDoseDat << '\n';
		kfacsByDoseDat << '\n';
		
		std::string starOut = outDir+ZIO::itoa(t)+".star";
		
		mdt.write(starOut);
		tomogramSetOut.setDose(t, starOut);*/
	}
	
	tomogramSetOut.write(outDir + "tomograms.star");
}
