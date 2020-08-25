
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/membrane/blob_3d.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/membrane/quadric_disc_fit.h>
#include <src/jaz/membrane/phaseline_average.h>

#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{	
	IOParser parser;
	
	std::string outTag, tomoListFn, catFnN, catFnS, catFnC;
	
	int num_threads;
    double diameter, binning, thickness;
	bool diag;
		
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		catFnN = parser.getOption("--co", "Catalogue file, outer membrane particles");
		catFnS = parser.getOption("--ci", "Catalogue file, inner membrane particles");
		catFnC = parser.getOption("--cc", "Catalogue file, substack particles");
		
		tomoListFn = parser.getOption("--t", "Tomogram list", "tomolist.txt");
		
		diameter = textToDouble(parser.getOption("--d", "Disc diameter in binned pixels", "64"));
		thickness = textToDouble(parser.getOption("--th", "Membrane thickness in binned pixels", "3"));
		
				
		binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
		diag = parser.checkOption("--diag", "Write out diagnostic information");
		
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		outTag = parser.getOption("--o", "Output filename pattern");
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	Catalogue catN(catFnN), catS(catFnS), catC(catFnC);			
	TomoList tomoList(tomoListFn);
	
	const int tc = catN.particles.size();
	
	if (catN.particles.size() != catS.particles.size())
	{
		REPORT_ERROR_STR("Unequal numbers of particles in "
						 << catFnN << " and " << catFnS << ".");
	}
		
	int dummy = std::system(("mkdir -p " + outTag).c_str());
	
	std::string diagTag = outTag + "/diag/";
	
	if (diag)
	{
		int dummy = std::system(("mkdir -p " + diagTag).c_str());
	}
	
	for (int tt = 0; tt < tc; tt++)
	{
		DynamoParticle& partC = catC.particles[tt];
		
		DynamoParticle& partN = catN.particles[tt];
		DynamoParticle& partS = catS.particles[tt];
		
		const int t = partC.tomo;
		
		if (!tomoList.isKnown(t))
		{
			REPORT_ERROR_STR("Bad tomogram: " << t << ".");
		}
						
		std::string stackFn = tomoList.getTiltSeriesFilename(t);
		std::string projFn = tomoList.getProjectionsFilename(t);
				
		BufferedImage<float> stack;
		stack.read(stackFn);
		
		const int  w = stack.xdim;
		const int  h = stack.ydim;
		const int fc = stack.zdim;
		
		int w0, h0, d0;
		std::vector<d4Matrix> projTomo = ProjectionIO::read(projFn, w0, h0, d0);
	
		QuadricDiscFit mdf(
			stack, projTomo, 
			partN.getPosition(),
			partS.getPosition(),
			binning, diameter, thickness, num_threads);
		
		std::vector<double> x0 = mdf.getInitial();
		
		std::vector<double> x1 = NelderMead::optimize(x0, mdf, 1e-5, 1e-6, 10000);
		
		std::cout << "final: " << mdf.f(x1, 0) << std::endl;
		
		BufferedImage<float> phaselines = mdf.writePhaselines(x1, w, h, projTomo);
		
		BufferedImage<float> stackAvg(w,h,fc);
		
		for (int f = 0; f < fc; f++)
		{
			std::vector<double> avg = PhaseLineAverage::averageNN(
				stack.getSliceRef(f), phaselines.getSliceRef(f), w);
			
			BufferedImage<float> expa = PhaseLineAverage::expandLIN(
				stack.getSliceRef(f), phaselines.getSliceRef(f), avg);
			
			stackAvg.copySliceFrom(f, expa);
		}
			
		stack -= stackAvg;
		
		if (diag)
		{
			std::string tag = diagTag + partC.getFormattedTag();
			
			mdf.localTomo.write(tag + "_3D_tomo.mrc");
					
			BufferedImage<float> vis0 = mdf.visualize(x0);
			BufferedImage<float> vis1 = mdf.visualize(x1);
			
			vis0.write(tag + "_3D_initial.mrc");
			vis1.write(tag + "_3D_final.mrc");
			
			phaselines.write(tag + "_phase_lines.mrc");
			
			BufferedImage<float> vis2D = mdf.visualize2D(x1, stack, projTomo);
			vis2D.write(tag + "_fit.mrc");
			
			stackAvg.write(tag + "_average.mrc");
		}
		
		stack.write(outTag + "/" + partC.getFormattedTag() + ".mrc");
	}
	
	return 0;
}

