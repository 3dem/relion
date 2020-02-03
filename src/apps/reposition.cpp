#include <src/jaz/obs_model.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/reference_map.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/interpolation.h>
#include <src/jaz/new_ft.h>
#include <src/ctf.h>
#include <src/metadata_table.h>
#include <omp.h>


using namespace gravis;

int main(int argc, char *argv[])
{
	ReferenceMap reference;
	
	IOParser parser;
	
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	
	std::string starFn = parser.getOption("--i", "Input STAR file containing the particles");
	const int max_shift = textToInteger(parser.getOption("--max_shift", "Maximal allowed shift", "10"));
	const double pad = textToDouble(parser.getOption("--cc_pad", "Cross-correlation padding", "1"));
	
	reference.read(parser, argc, argv);
	
	const int nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
	
	std::string outFn = parser.getOption("--o", "Output path", "repositioned/");
	
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	reference.load(1, true);
	reference.projectors[0].interpolator = TRILINEAR;
	
	// Get dimensions
	const int s = reference.s;
	const int sh = s/2 + 1;
	
	ObservationModel obsModel;	
	MetaDataTable mdt0;
	
	ObservationModel::loadSafely(starFn, obsModel, mdt0);
	
	std::vector<MetaDataTable> allMdts = StackHelper::splitByMicrographName(mdt0);
	
	const int mgc = allMdts.size();
	
	
	Image<RFLOAT> freqWgh(sh,s);
	freqWgh.data.initConstant(1.0);
	
	freqWgh = FilterHelper::raisedCosEnvFreq2D(
				freqWgh, reference.k_out-1, reference.k_out+1);
	
	std::vector<Image<dComplex>> ccsFs(nr_omp_threads);
	std::vector<Image<double>> ccsRs(nr_omp_threads);
	
	const int s2 = (int)(pad * s);
	const int sh2 = s2/2 + 1;
	
	for (int t = 0; t < nr_omp_threads; t++)
	{
		ccsFs[t] = Image<dComplex>(sh2,s2);
		ccsRs[t] = Image<double>(s2,s2);
	}
	
	NewFFT::DoublePlan plan(s2,s2,1);
	
	MetaDataTable outMdt;
	
	for (int m = 0; m < mgc; m++)
	{
		std::cout << m << " / " << mgc << "\n";
		
		const int pc = allMdts[m].numberOfObjects();
		
		std::vector<Image<Complex>> pred = reference.predictAll(
					allMdts[m], obsModel, ReferenceMap::Own, nr_omp_threads,
					true, true, false);
		
		std::vector<Image<Complex>> obs = StackHelper::loadStackFS(
					allMdts[m], "", nr_omp_threads, true, &obsModel);
		
		#pragma omp parallel for num_threads(nr_omp_threads)
		for (int p = 0; p < pc; p++)
		{
			int t = omp_get_thread_num();
			
			int opticsGroup = obsModel.getOpticsGroup(allMdts[m], p);		
			double angpix = obsModel.getPixelSize(opticsGroup);
			
			for (int y = 0; y < s; y++)
			for (int x = 0; x < sh; x++)
			{
				Complex z = obs[p](y,x) * freqWgh(y,x) * pred[p](y,x).conj();
				
				const int yy = y < sh? y : s2 - (s - y);
				
				ccsFs[t](yy,x) = dComplex(z.real, z.imag);
			}
			
			NewFFT::inverseFourierTransform(ccsFs[t](), ccsRs[t](), plan);
			
			d2Vector shift = Interpolation::quadraticMaxWrapXY(
						ccsRs[t], 1e-12, pad * max_shift, pad * max_shift);
			
			if (shift.x >= sh2) shift.x -= s2;
			if (shift.y >= sh2) shift.y -= s2;
						
			double xoff, yoff;
			
			allMdts[m].getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, p);
			allMdts[m].getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, p);
			
			xoff /= angpix;
			yoff /= angpix;
			
			xoff -= shift.x;
			yoff -= shift.y;
			
			xoff *= angpix;
			yoff *= angpix;
			
			allMdts[m].setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, p);
			allMdts[m].setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, p);
		}
		
		outMdt.append(allMdts[m]);
	}
	
	if (outFn.length() > 0 && outFn[outFn.length()-1] != '/')
	{
		outFn = outFn + "/";
	}
	
	std::string command = " mkdir -p " + outFn;
	int res = system(command.c_str());
	
	outMdt.write(outFn+"repositioned.star");
	
	return RELION_EXIT_SUCCESS;
}
