#include "magnification_estimator.h"
#include "magnification_helper.h"
#include "ctf_refiner.h"
#include "equation2x2.h"

#include <src/jaz/reference_map.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/vtk_helper.h>

#include <omp.h>


using namespace gravis;

MagnificationEstimator::MagnificationEstimator()
{
	
}

void MagnificationEstimator::read(
		IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_mag", 
				"Inner freq. threshold for anisotropic magnification estimation [Angst]", "20.0"));
}

void MagnificationEstimator::init(
		int verb, int s, int nr_omp_threads, 
		bool debug, bool diag, 
		std::string outPath, 
		ReferenceMap* reference, 
		ObservationModel* obsModel)
{
	this->verb = verb;
	this->s = s;
	sh = s/2 + 1;
	this->nr_omp_threads = nr_omp_threads;
	
	this->debug = debug;
	this->diag = diag;
	this->outPath = outPath;
	
	this->reference = reference;
	this->obsModel = obsModel;
	
	angpix = obsModel->getPixelSize(0);
	
	ready = true;
}

void MagnificationEstimator::processMicrograph(
		long g, MetaDataTable& mdt, 
		const std::vector<Image<Complex>>& obs, 
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR_STR("ERROR: MagnificationEstimator::processMicrograph: "
						 << "MagnificationEstimator not initialized.");
	}
	
	std::vector<Volume<Equation2x2>> magEqs(nr_omp_threads);
	
	for (int i = 0; i < nr_omp_threads; i++)
	{
		magEqs[i] = Volume<Equation2x2>(sh,s,1);
	}
	
	const int pc = mdt.numberOfObjects();
	
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		CTF ctf;
		ctf.readByGroup(mdt, obsModel->opticsMdt, p);
		
		int threadnum = omp_get_thread_num();
	
		MagnificationHelper::updateScaleFreq(pred[p], obs[p], ctf, angpix, magEqs[threadnum]);
	}
	
	for (int i = 1; i < nr_omp_threads; i++)
	{
		magEqs[0] += magEqs[i];
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
		
	MagnificationHelper::writeEQs(magEqs[0], outRoot+"_mag");
}

void MagnificationEstimator::parametricFit(
		const std::vector<MetaDataTable> &mdts, MetaDataTable &mdtOut)
{
	if (!ready)
	{
		REPORT_ERROR_STR("ERROR: MagnificationEstimator::parametricFit: "
					 << "MagnificationEstimator not initialized.");
	}
	
	if (verb > 0)
	{
		std::cout << " + Fitting anisotropic magnification ..." << std::endl;
	}
	
	Volume<Equation2x2> magEqs(sh,s,1), magEqsG(sh,s,1);
	
	const int gc = mdts.size();
				
	for (long g = 0; g < gc; g++)
	{
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
		
		MagnificationHelper::readEQs(outRoot+"_mag", magEqsG);
		
		magEqs += magEqsG;
	}
	
	Image<RFLOAT> flowx, flowy;
	MagnificationHelper::solvePerPixel(magEqs, flowx, flowy);
	
	if (debug)
	{
		flowx.write(outPath+"_mag_disp_x.mrc");
		flowy.write(outPath+"_mag_disp_y.mrc");
	}
	
	Image<RFLOAT> flowxFull, flowyFull;
	FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
	FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);
	
	if (debug)
	{
		VtkHelper::writeVTK(flowxFull, outPath+"_mag_disp_x.vtk");
		VtkHelper::writeVTK(flowyFull, outPath+"_mag_disp_y.vtk");
	}
	
	Image<RFLOAT> freqWght = reference->getHollowWeight(kmin);
	
	Image<RFLOAT> mat;
	MagnificationHelper::solveLinearlyFreq(magEqs, freqWght, mat, flowx, flowy);
	
	if (debug)
	{
		flowx.write(outPath+"_mag_disp_x_fit.mrc");
		flowy.write(outPath+"_mag_disp_y_fit.mrc");
		mat.write(outPath+"_mag_disp_mat.mrc");
	
		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);
		
		VtkHelper::writeVTK(flowxFull, outPath+"_mag_disp_x_fit.vtk");
		VtkHelper::writeVTK(flowyFull, outPath+"_mag_disp_y_fit.vtk");
	}
	
	std::ofstream os(outPath+"_mag_mat.txt");
	os << DIRECT_A2D_ELEM(mat.data, 0, 0) << " " << DIRECT_A2D_ELEM(mat.data, 0, 1) << "\n";
	os << DIRECT_A2D_ELEM(mat.data, 1, 0) << " " << DIRECT_A2D_ELEM(mat.data, 1, 1) << "\n";
	os.close();
}

bool MagnificationEstimator::isFinished(
		const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::isFinished: DefocusEstimator not initialized.");
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	return exists(outRoot+"_mag_Axx.mrc")
	    && exists(outRoot+"_mag_Axy.mrc")
		&& exists(outRoot+"_mag_Ayy.mrc")
		&& exists(outRoot+"_mag_bx.mrc")
		&& exists(outRoot+"_mag_by.mrc");
}
