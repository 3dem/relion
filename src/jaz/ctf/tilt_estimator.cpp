#include "tilt_estimator.h"
#include "tilt_helper.h"
#include "ctf_refiner.h"

#include <src/jaz/obs_model.h>
#include <src/jaz/reference_map.h>
#include <src/jaz/image_op.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/gravis/t2Vector.h>

#include <src/args.h>
#include <src/ctf.h>

#include <set>

#include <omp.h>

using namespace gravis;


TiltEstimator::TiltEstimator()
	:	ready(false)
{}

void TiltEstimator::read(IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_tilt", 
		"Inner freq. threshold for beamtilt estimation [Angst]", "20.0"));
	
	std::string aberrToken = "--odd_aberr_max_n";
	
	aberr_n_max  = textToInteger(parser.getOption(aberrToken, 
		"Maximum degree of Zernike polynomials used to fit odd (i.e. antisymmetrical) aberrations", "0"));
}

void TiltEstimator::init(
		int verb, int s, int nr_omp_threads, 
		bool debug, bool diag, std::string outPath,
		ReferenceMap* reference, ObservationModel* obsModel)
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

void TiltEstimator::processMicrograph(
		long g, MetaDataTable& mdt, 
		const std::vector<Image<Complex>>& obs,
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::processMicrograph: TiltEstimator not initialized.");
	}
	
	const int pc = mdt.numberOfObjects();
	
	std::vector<int> optGroups = obsModel->getOptGroupsPresent(mdt);	
	const int cc = optGroups.size();
	
	std::vector<int> groupToIndex(obsModel->numberOfOpticsGroups()+1, -1);
	
	for (int i = 0; i < cc; i++)
	{
		groupToIndex[optGroups[i]] = i;
	}
		
	std::vector<Image<Complex>> xyAcc(nr_omp_threads*cc);
	std::vector<Image<RFLOAT>> wAcc(nr_omp_threads*cc);
	
	for (int i = 0; i < nr_omp_threads*cc; i++)
	{
		xyAcc[i] = Image<Complex>(sh,s);
		xyAcc[i].data.initZeros();
		
		wAcc[i] = Image<RFLOAT>(sh,s);
		wAcc[i].data.initZeros();
	}
	
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		CTF ctf;
		ctf.readByGroup(mdt, obsModel, p);
		
		int threadnum = omp_get_thread_num();
		
		int og;
		mdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);

		const int ci = groupToIndex[og];

		TiltHelper::updateTiltShift(
			pred[p], obs[p], ctf, angpix, 
			xyAcc[cc*threadnum + ci], 
			wAcc[cc*threadnum + ci]);
	}
	
	// Combine the accumulated weights from all threads for this subset, 
	// store weighted sums in xyAccSum and wAccSum
	
	for (int ci = 0; ci < cc; ci++)
	{
		Image<Complex> xyAccSum(sh,s);
		Image<RFLOAT> wAccSum(sh,s);
		
		for (int threadnum = 0; threadnum < nr_omp_threads; threadnum++)
		{
			ImageOp::linearCombination(xyAccSum, xyAcc[cc*threadnum + ci], 1.0, 1.0, xyAccSum);
			ImageOp::linearCombination(wAccSum, wAcc[cc*threadnum + ci], 1.0, 1.0, wAccSum);
		}
		
		// Write out the intermediate results per-micrograph:
		
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
		
		std::stringstream sts;
		sts << optGroups[ci];
		
		ComplexIO::write(xyAccSum(), outRoot + "_xyAcc_optics-class_" + sts.str(), ".mrc");
		wAccSum.write(outRoot+"_wAcc_optics-class_" + sts.str() + ".mrc");
	}
}

void TiltEstimator::parametricFit(
		const std::vector<MetaDataTable>& mdts, 
		MetaDataTable& optOut)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::parametricFit: TiltEstimator not initialized.");
	}
	
	if (verb > 0)
	{
		std::cout << " + Fitting beam tilt ..." << std::endl;
	}
	
	const int gc = mdts.size();	
	const int ogc = obsModel->numberOfOpticsGroups();
	
	std::vector<d2Vector> tiltPerClass(ogc, d2Vector(0.0, 0.0));
	
	std::vector<bool> groupUsed(ogc,false);
	
	for (int og = 0; og < ogc; og++)
	{	
		double Cs = obsModel->Cs[og];
		double lambda = obsModel->lambda[og];
		
		std::stringstream sts;
		sts << og+1;
		std::string cns = sts.str();
		
		Image<Complex> xyAccSum(sh,s);
		Image<RFLOAT> wAccSum(sh,s);
		
		xyAccSum.data.initZeros();
		wAccSum.data.initZeros();
			
		for (long g = 0; g < gc; g++)
		{
			std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
			
			if (   exists(outRoot+"_xyAcc_optics-class_"+cns+"_real.mrc")
				&& exists(outRoot+"_xyAcc_optics-class_"+cns+"_imag.mrc")
				&& exists(outRoot+ "_wAcc_optics-class_"+cns+".mrc"))
			{
				Image<Complex> xyAcc;
				Image<RFLOAT> wAcc;
				
				wAcc.read(outRoot+"_wAcc_optics-class_"+cns+".mrc");
				ComplexIO::read(xyAcc, outRoot+"_xyAcc_optics-class_"+cns, ".mrc");
				
				xyAccSum() += xyAcc();
				wAccSum()  +=  wAcc();
				
				groupUsed[og] = true;
			}	
		}
		
		if (!groupUsed[og])
		{
			continue;
		}
			
		Image<RFLOAT> wgh, phase, fit, phaseFull, fitFull;
		
		FilterHelper::getPhase(xyAccSum, phase);
		
		Image<Complex> xyNrm(sh,s);
	
		double kmin_px = obsModel->angToPix(kmin, s, og);
				 
		Image<RFLOAT> wgh0 = reference->getHollowWeight(kmin_px);
	
		FilterHelper::multiply(wAccSum, wgh0, wgh);
		
		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			xyNrm(y,x) = wAccSum(y,x) > 0.0? xyAccSum(y,x)/wAccSum(y,x) : Complex(0.0, 0.0);
		}
		
		if (debug)
		{
			Image<RFLOAT> wghFull;
			FftwHelper::decenterDouble2D(wgh(), wghFull());
			
			ImageLog::write(wghFull, outPath + "beamtilt_weight-full_optics-class_"+cns);
		}
		
		FftwHelper::decenterUnflip2D(phase.data, phaseFull.data);
			
		ImageLog::write(phaseFull, outPath + "beamtilt_delta-phase_per-pixel_optics-class_"+cns);
		
		double shift_x(0), shift_y(0), tilt_x(0), tilt_y(0);
		
		if (aberr_n_max < 3)
		{
			TiltHelper::fitTiltShift(
				phase, wgh, Cs, lambda, angpix,
				&shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
				
			FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
			
			ImageLog::write(fitFull, outPath + "beamtilt_delta-phase_lin-fit_optics-class_"+cns);
						
			TiltHelper::optimizeTilt(
					xyNrm, wgh, Cs, lambda, angpix, false,
					shift_x, shift_y, tilt_x, tilt_y,
					&shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
			
			FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
			
			ImageLog::write(fitFull, outPath+"beamtilt_delta-phase_iter-fit_optics-class_"+cns);
			
			optOut.setValue(EMDL_IMAGE_BEAMTILT_X, tilt_x, og);
			optOut.setValue(EMDL_IMAGE_BEAMTILT_Y, tilt_y, og);
		}
		else
		{
			Image<RFLOAT> one(sh,s);
			one.data.initConstant(1);
			
			std::vector<double> Zernike_coeffs = TiltHelper::fitOddZernike(
						xyNrm, wgh, angpix, aberr_n_max, &fit);
						
			FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
			
			std::stringstream sts;
			sts << aberr_n_max;
			
			ImageLog::write(fitFull, outPath + "beamtilt_delta-phase_lin-fit_optics-class_"
							+cns+"_N-"+sts.str());
			
			if (debug)
			{
				Image<RFLOAT> residual;
				residual.data = phaseFull.data - fitFull.data;
				
				ImageLog::write(residual, outPath + "beamtilt_delta-phase_lin-fit_optics-class_"
								+cns+"_N-"+sts.str()+"_residual");
			}
			
			Image<RFLOAT> plot0 = TiltHelper::plotOddZernike(Zernike_coeffs, s, angpix);
			ImageLog::write(plot0, outPath + "debug_Z0");
			
			TiltHelper::extractTilt(Zernike_coeffs, tilt_x, tilt_y, Cs, lambda);
			
			/*
			{
				Image<RFLOAT> plot1 = TiltHelper::plotOddZernike(Zernike_coeffs, s, angpix);
				ImageLog::write(plot1, outPath + "debug_Z1");
				
				Image<RFLOAT> plot2 = TiltHelper::plotTilt(tilt_x, tilt_y, s, angpix, Cs, lambda);
				ImageLog::write(plot2, outPath + "debug_t");
				
				Image<RFLOAT> plot3;
				plot3() = plot1() - plot2();
				ImageLog::write(plot3, outPath + "debug_t+Z1");
			}
			*/
			
			optOut.setValue(EMDL_IMAGE_BEAMTILT_X, tilt_x, og);
			optOut.setValue(EMDL_IMAGE_BEAMTILT_Y, tilt_y, og);
			optOut.setValue(EMDL_IMAGE_ODD_ZERNIKE_COEFFS, Zernike_coeffs, og);
		}
	}
}

bool TiltEstimator::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::isFinished: DefocusEstimator not initialized.");
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	bool allDone = true;
	
	const int ogc = obsModel->numberOfOpticsGroups();
	
	for (int og = 0; og < ogc; og++)
	{	
		std::stringstream sts;
		sts << og+1;
		std::string ogs = sts.str();
		
		if (   !exists(outRoot+"_xyAcc_optics-class_"+ogs+"_real.mrc")
			|| !exists(outRoot+"_xyAcc_optics-class_"+ogs+"_imag.mrc")
			|| !exists(outRoot+"_wAcc_optics-class_"+ogs+".mrc"))
		{
			allDone = false;
			break;
		}
	}
	
	return allDone;
}
