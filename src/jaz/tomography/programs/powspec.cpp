#include "powspec.h"
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/membrane/blob_3d.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/conversion.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>
#include <src/jaz/util/image_file_helper.h>

#include <omp.h>

using namespace gravis;


void PowspecProgram::run()
{
	std::istringstream sts(stackFn);
	std::string token;
	
	std::vector<std::string> fns;
	
	while(std::getline(sts, token, ','))
	{
		fns.push_back(token);
	}
	
	const int outSize = res/2 + 1;
	
	std::vector<double> out;
	
	
	std::vector<double> ps0;
	
	const bool do_ratio = ratio >= 0;
	
	if (!separate)
	{
		out = std::vector<double>(outSize, 0.0);
	}
			
	for (int i = 0; i < fns.size(); i++)
	{
		std::cout << "processing " << fns[i] << "..." << std::endl;
		
		BufferedImage<float> ts;
		ts.read(fns[i]);
		
		if (do_ratio)
		{
			ps0 = PowerSpectrum::periodogramAverage1D(ts, res, 2.0, ratio);
		}
		
		const double scale = (double) (res * res) / (double) fns.size() * (double) ts.zdim;
		
		for (int f = 0; f < ts.zdim; f++)
		{	
			std::cout << "    frame " << f << std::endl;
			
			if (separate)
			{
				out = std::vector<double>(outSize, 0.0);
			}
			
			std::vector<double> ps = PowerSpectrum::periodogramAverage1D(ts, res, 2.0, f);
			
			for (int i = 0; i < outSize; i++)
			{
				out[i] += scale * ps[i];
			}
			
			if (separate)
			{
				std::stringstream sts;
				sts << f;
				
				if (do_ratio)
				{			
					std::string outFnR = outFn + "_frame_" + sts.str() + "_ratio.dat";				
					std::ofstream ofsR(outFnR);
					
					for (int i = 1; i < outSize; i++)
					{
						ofsR << i / (double)(res * pixSize)  << " " << ps[i] / ps0[i] << "\n";
					}
				}
				
				std::string outFnF = outFn + "_frame_" + sts.str() + ".dat";				
				std::ofstream ofs(outFnF);
				
				for (int i = 1; i < outSize; i++)
				{
					ofs << i / (double)(res * pixSize)  << " " << out[i] << "\n";
				}
			}
		}
		
		std::cout << std::endl;
	}
	
	if (!separate)
	{
		std::ofstream ofs(outFn);
		
		for (int i = 1; i < outSize; i++)
		{
			ofs << i / (double)(res * pixSize)  << " " << out[i] << "\n";
		}
	}
	
}
