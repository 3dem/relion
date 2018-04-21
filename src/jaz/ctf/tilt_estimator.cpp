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

#include <src/args.h>
#include <src/ctf.h>

#include <omp.h>

TiltEstimator::TiltEstimator()
	:	ready(false)
{}

void TiltEstimator::read(IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_tilt", 
						"Inner freq. threshold for beamtilt estimation [Angst]", "20.0"));
		
	aniso = parser.checkOption("--anisotropic_tilt", "Use anisotropic coma model");
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
	
	angpix = obsModel->angpix;
	
	ready = true;
}

void TiltEstimator::processMicrograph(
		long g, MetaDataTable& mdt, 
		const std::vector<Image<Complex>>& obs)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::processMicrograph: TiltEstimator not initialized.");
	}
	
	const int pc = mdt.numberOfObjects();
	
	std::vector<Image<Complex>> pred = reference->predictAll(
		mdt, *obsModel, ReferenceMap::Opposite, nr_omp_threads, false, false);
	
	std::vector<Image<Complex>> xyAcc(nr_omp_threads);
	std::vector<Image<RFLOAT>> wAcc(nr_omp_threads);
	
	for (int i = 0; i < nr_omp_threads; i++)
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
		ctf.read(mdt, mdt, p);
		
		int threadnum = omp_get_thread_num();
		
		TiltHelper::updateTiltShift(
			pred[p], obs[p], ctf, obsModel->angpix, 
			xyAcc[threadnum], wAcc[threadnum]);
	}
	
	// Combine the accumulated weights from all threads for this subset, 
	// store weighted sums in xyAccSum and wAccSum
	
	Image<Complex> xyAccSum(sh,s);
	Image<RFLOAT> wAccSum(sh,s);
	
	for (int i = 0; i < nr_omp_threads; i++)
	{
		ImageOp::linearCombination(xyAccSum, xyAcc[i], 1.0, 1.0, xyAccSum);
		ImageOp::linearCombination(wAccSum, wAcc[i], 1.0, 1.0, wAccSum);
	}
	
	// Write out the intermediate results per-micrograph:
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	ComplexIO::write(xyAccSum(), outRoot+"_xyAcc", ".mrc");
	wAccSum.write(outRoot+"_wAcc.mrc");
}

void TiltEstimator::parametricFit(
		const std::vector<MetaDataTable>& mdts, 
		double Cs, double lambda, 
		MetaDataTable& mdtOut)
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
		
	Image<Complex> xyAccSum(sh,s);
	Image<RFLOAT> wAccSum(sh,s);
	
	xyAccSum.data.initZeros();
	wAccSum.data.initZeros();
		
	for (long g = 0; g < gc; g++)
	{
		Image<Complex> xyAcc;
		Image<RFLOAT> wAcc;
		
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
		
		wAcc.read(outRoot+"_wAcc.mrc");
		ComplexIO::read(xyAcc, outRoot+"_xyAcc", ".mrc");
		
		xyAccSum() += xyAcc();
		wAccSum() += wAcc();
	}
		
	Image<RFLOAT> wgh, phase, fit, phaseFull, fitFull;
	
	FilterHelper::getPhase(xyAccSum, phase);
	
	Image<Complex> xyNrm(sh,s);

	double kmin_px = obsModel->angToPix(kmin, s);
			 
	Image<RFLOAT> wgh0 = reference->getHollowWeight(kmin_px);

	FilterHelper::multiply(wAccSum, wgh0, wgh);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		xyNrm(y,x) = wAccSum(y,x) > 0.0? xyAccSum(y,x)/wAccSum(y,x) : Complex(0.0, 0.0);
	}
	
	Image<RFLOAT> wghFull;
	FftwHelper::decenterDouble2D(wgh(), wghFull());
	
	if (debug)
	{
		ImageLog::write(wghFull, outPath + "beamtilt_weight-full");
	}
	
	FftwHelper::decenterUnflip2D(phase.data, phaseFull.data);
		
	ImageLog::write(phaseFull, outPath + "beamtilt_delta-phase_per-pixel");
	
	double shift_x, shift_y, tilt_x, tilt_y;
	
	TiltHelper::fitTiltShift(
		phase, wgh, Cs, lambda, angpix,
		&shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
		
	FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
	
	ImageLog::write(fitFull, outPath + "beamtilt_delta-phase_lin-fit");
	
	std::ofstream os(outPath+"beamtilt_0.txt");
	os << "beamtilt_x = " << tilt_x << "\n";
	os << "beamtilt_y = " << tilt_y << "\n";
	os.close();
	
	double tilt_xx, tilt_xy, tilt_yy;
	
	if (aniso)
	{
		TiltHelper::optimizeAnisoTilt(
			xyNrm, wgh, Cs, lambda, angpix, false,
			shift_x, shift_y, tilt_x, tilt_y,
			&shift_x, &shift_y, &tilt_x, &tilt_y,
			&tilt_xx, &tilt_xy, &tilt_yy, &fit);
	}
	else
	{
		TiltHelper::optimizeTilt(
			xyNrm, wgh, Cs, lambda, angpix, false,
			shift_x, shift_y, tilt_x, tilt_y,
			&shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
	}
	
	FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
	
	ImageLog::write(fitFull, outPath+"beamtilt_delta-phase_iter-fit");
	
	std::ofstream os2(outPath+"fit_beamtilt_1.txt");
	os2 << "beamtilt_x = " << tilt_x << "\n";
	os2 << "beamtilt_y = " << tilt_y << "\n";
	os2.close();
		
	// Now write the beamtilt into mdtOut	
	const long tpc = mdtOut.numberOfObjects();
	
	for (long p = 0; p < tpc; p++)
	{
		mdtOut.setValue(EMDL_IMAGE_BEAMTILT_X, tilt_x, p);
		mdtOut.setValue(EMDL_IMAGE_BEAMTILT_Y, tilt_y, p);
	}
}

bool TiltEstimator::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::isFinished: DefocusEstimator not initialized.");
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	return exists(outRoot+"_xyAcc_real.mrc")
	    && exists(outRoot+"_xyAcc_imag.mrc")
	    && exists(outRoot+"_wAcc.mrc");
}
