#include "aberration_estimator.h"
#include "tilt_helper.h"
#include "ctf_refiner.h"

#include <src/jaz/obs_model.h>
#include <src/jaz/reference_map.h>
#include <src/jaz/img_proc/image_op.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/gravis/t2Vector.h>

#include <src/args.h>
#include <src/ctf.h>

#include <set>

#include <omp.h>

using namespace gravis;


AberrationEstimator::AberrationEstimator()
	:	ready(false)
{}

void AberrationEstimator::read(IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_aberr", 
		"Inner freq. threshold for symmetrical aberration estimation [Angst]", "20.0"));
	
	std::string aberrToken = "--even_aberr_max_n";
	
	aberr_n_max  = textToInteger(parser.getOption(aberrToken, 
		"Maximum degree of Zernike polynomials used to fit even (i.e. symmetrical) aberrations", "4"));
	
	xring0 = textToDouble(parser.getOption("--xr0_a", "Exclusion ring start (A)", "-1"));	
	xring1 = textToDouble(parser.getOption("--xr1_a", "Exclusion ring end (A)", "-1"));
}

void AberrationEstimator::init(
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

void AberrationEstimator::processMicrograph(
		long g, MetaDataTable& mdt, 
		const std::vector<Image<Complex>>& obs,
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: AberrationEstimator::processMicrograph: AberrationEstimator not initialized.");
	}
	
	const int pc = mdt.numberOfObjects();
	
	std::vector<int> optGroups = obsModel->getOptGroupsPresent(mdt);	
	const int cc = optGroups.size();
	
	std::vector<int> groupToIndex(obsModel->numberOfOpticsGroups()+1, -1);
	
	for (int i = 0; i < cc; i++)
	{
		groupToIndex[optGroups[i]] = i;
	}
	
	std::vector<Image<RFLOAT>>
		Axx(nr_omp_threads*cc, Image<RFLOAT>(sh,s)),
		Axy(nr_omp_threads*cc, Image<RFLOAT>(sh,s)),
		Ayy(nr_omp_threads*cc, Image<RFLOAT>(sh,s)),
		bx(nr_omp_threads*cc, Image<RFLOAT>(sh,s)),
		by(nr_omp_threads*cc, Image<RFLOAT>(sh,s));
	
	const double as = (double)s * angpix;
	
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		CTF ctf;
		ctf.readByGroup(mdt, obsModel, p);
		
		int threadnum = omp_get_thread_num();
		
		int og;
		mdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);

		const int ci = groupToIndex[og];
		
		const int t = cc * threadnum + ci;
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double xf = x;
			const double yf = y < sh? y : y - s;
			
			const double gamma_i = ctf.getGamma(xf/as, yf/as);
			const double cg = cos(gamma_i);
			const double sg = sin(gamma_i);

			Complex zobs = obs[p](y,x);
			Complex zprd = pred[p](y,x);

			double zz = zobs.real*zprd.real + zobs.imag*zprd.imag;
			double nr = zprd.norm();

			Axx[t](y,x) += nr * sg * sg;
			Axy[t](y,x) += nr * cg * sg;
			Ayy[t](y,x) += nr * cg * cg;

			bx[t](y,x) -= zz * sg;
			by[t](y,x) -= zz * cg;
		}
	}
	
	// Combine the accumulated weights from all threads for this subset
	
	for (int ci = 0; ci < cc; ci++)
	{
		Image<RFLOAT> 
			AxxSum(sh,s), AxySum(sh,s), AyySum(sh,s),
			bxSum(sh,s), bySum(sh,s);
		
		for (int threadnum = 0; threadnum < nr_omp_threads; threadnum++)
		{
			ImageOp::linearCombination(AxxSum, Axx[cc*threadnum + ci], 1.0, 1.0, AxxSum);
			ImageOp::linearCombination(AxySum, Axy[cc*threadnum + ci], 1.0, 1.0, AxySum);
			ImageOp::linearCombination(AyySum, Ayy[cc*threadnum + ci], 1.0, 1.0, AyySum);
			
			ImageOp::linearCombination(bxSum, bx[cc*threadnum + ci], 1.0, 1.0, bxSum);
			ImageOp::linearCombination(bySum, by[cc*threadnum + ci], 1.0, 1.0, bySum);
		}
		
		// Write out the intermediate results per-micrograph:
		
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
		
		std::stringstream sts;
		sts << optGroups[ci];
		
		AxxSum.write(outRoot+"_aberr-Axx_optics-class_" + sts.str() + ".mrc");
		AxySum.write(outRoot+"_aberr-Axy_optics-class_" + sts.str() + ".mrc");
		AyySum.write(outRoot+"_aberr-Ayy_optics-class_" + sts.str() + ".mrc");
		
		bxSum.write(outRoot+"_aberr-bx_optics-class_" + sts.str() + ".mrc");
		bySum.write(outRoot+"_aberr-by_optics-class_" + sts.str() + ".mrc");
	}
}

void AberrationEstimator::parametricFit(
		const std::vector<MetaDataTable>& mdts, 
		MetaDataTable& optOut)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: AberrationEstimator::parametricFit: AberrationEstimator not initialized.");
	}
	
	if (verb > 0)
	{
		std::cout << " + Fitting symmetrical aberrations ..." << std::endl;
	}
	
	const int gc = mdts.size();	
	const int ogc = obsModel->numberOfOpticsGroups();
		
	std::vector<bool> groupUsed(ogc,false);
	
	for (int og = 0; og < ogc; og++)
	{	
		std::stringstream sts;
		sts << og+1;
		std::string cns = sts.str();
		
		Image<RFLOAT> 
			AxxSum(sh,s), AxySum(sh,s), AyySum(sh,s),
			bxSum(sh,s), bySum(sh,s);
		
		for (long g = 0; g < gc; g++)
		{
			std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
			
			if (   exists(outRoot+"_aberr-Axx_optics-class_" + sts.str() + ".mrc")
				&& exists(outRoot+"_aberr-Axy_optics-class_" + sts.str() + ".mrc")
				&& exists(outRoot+"_aberr-Ayy_optics-class_" + sts.str() + ".mrc")
				&& exists(outRoot+"_aberr-bx_optics-class_" + sts.str() + ".mrc")
				&& exists(outRoot+"_aberr-by_optics-class_" + sts.str() + ".mrc"))
			{
				
				Image<RFLOAT> 
					Axx(sh,s), Axy(sh,s), Ayy(sh,s),
					bx(sh,s), by(sh,s);
				
				Axx.read(outRoot+"_aberr-Axx_optics-class_" + sts.str() + ".mrc");
				Axy.read(outRoot+"_aberr-Axy_optics-class_" + sts.str() + ".mrc");
				Ayy.read(outRoot+"_aberr-Ayy_optics-class_" + sts.str() + ".mrc");
				
				bx.read(outRoot+"_aberr-bx_optics-class_" + sts.str() + ".mrc");
				by.read(outRoot+"_aberr-by_optics-class_" + sts.str() + ".mrc");
				
				AxxSum() += Axx();
				AxySum() += Axy();
				AyySum() += Ayy();
				
				bxSum() += bx();
				bySum() += by();
				
				groupUsed[og] = true;
			}	
		}
		
		if (!groupUsed[og])
		{
			continue;
		}
			
		Image<RFLOAT> wgh(sh,s), phase(sh,s);
		Image<Complex> optXY(sh,s); 
				
		double kmin_px = obsModel->angToPix(kmin, s, og);
		wgh = reference->getHollowWeight(kmin_px);
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			d2Matrix A(
				AxxSum(y,x), AxySum(y,x),
				AxySum(y,x), AyySum(y,x));

			d2Vector b(bxSum(y,x), bySum(y,x));

			double det = A(0,0) * A(1,1) - A(1,0) * A(0,1);
			
			if (det != 0.0)
			{
				d2Matrix Ai = A;
				Ai.invert();

				d2Vector opt = Ai * b;

				optXY(y,x) = Complex(opt.x, opt.y);
				phase(y,x) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
				wgh(y,x) *= sqrt(sqrt(std::abs(det)));
			}
			else
			{
				optXY(y,x) = 0.0;
				phase(y,x) = 0.0;
				wgh(y,x) = 0.0;
			}
		}
		
		if (xring1 > 0.0)
		{
			for (int y = 0; y < s; y++)
			for (int x = 0; x < sh; x++)
			{
				double xx = x;
				double yy = y < sh? y : y - s;
				double rp = sqrt(xx*xx + yy*yy);
				double ra = s * angpix / rp;
				
				if (ra > xring0 && ra <= xring1)
				{
					wgh(y,x) = 0.0;
				}
			}
		}
				
		if (debug)
		{
			Image<RFLOAT> full;
			FftwHelper::decenterDouble2D(wgh(), full());
			ImageLog::write(full, outPath + "aberr_weight-full_optics-class_"+cns);
		}
		
		Image<RFLOAT> fit, phaseFull, fitFull;		
		FftwHelper::decenterDouble2D(phase.data, phaseFull.data);
		ImageLog::write(phaseFull, outPath + "aberr_delta-phase_per-pixel_optics-class_"+cns);
		
		
		
		
		{
			std::vector<double> Zernike_coeffs = TiltHelper::fitEvenZernike(
						phase, wgh, angpix, aberr_n_max, &fit);
						
			FftwHelper::decenterDouble2D(fit.data, fitFull.data);
			
			std::stringstream sts;
			sts << aberr_n_max;
			
			ImageLog::write(fitFull, outPath + "aberr_delta-phase_lin-fit_optics-class_"
							+cns+"_N-"+sts.str());
			if (debug)
			{
				Image<RFLOAT> residual;
				residual.data = phaseFull.data - fitFull.data;
				
				ImageLog::write(residual, outPath + "aberr_delta-phase_lin-fit_optics-class_"
								+cns+"_N-"+sts.str()+"_residual");
			}
			
			std::vector<double> Zernike_coeffs_opt = TiltHelper::optimiseEvenZernike(
						optXY, wgh, angpix, aberr_n_max, Zernike_coeffs, &fit);
				
			FftwHelper::decenterDouble2D(fit.data, fitFull.data);
						
			ImageLog::write(fitFull, outPath + "aberr_delta-phase_iter-fit_optics-class_"
							+cns+"_N-"+sts.str());
			
			// extract Q0, Cs, defocus and astigmatism?
			
			optOut.setValue(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, Zernike_coeffs_opt, og);
		}
	}
}

bool AberrationEstimator::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: AberrationEstimator::isFinished: AberrationEstimator not initialized.");
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	bool allDone = true;
	
	std::vector<int> ogs = obsModel->getOptGroupsPresent(mdt);
	
	for (int i = 0; i < ogs.size(); i++)
	{	
		const int og = ogs[i];
		
		std::stringstream sts;
		sts << og;
		
		if (   !exists(outRoot+"_aberr-Axx_optics-class_" + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-Axy_optics-class_" + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-Ayy_optics-class_" + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-bx_optics-class_"  + sts.str() + ".mrc")
			|| !exists(outRoot+"_aberr-by_optics-class_"  + sts.str() + ".mrc"))
		{
			allDone = false;
			break;
		}
	}
	
	return allDone;
}
